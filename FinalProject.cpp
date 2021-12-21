#include <iostream>
#include <fstream>
#include <math.h>
#include "matplotlibcpp.h"
#include <vector>
#include <queue>
#include <list>
#include <tuple>
using namespace std;
namespace plt = matplotlibcpp;

// functions explained below
tuple<map <int, vector<int>>, map <int, vector<int>>, vector<string>> iterCol (int node, int corrRow, int corrCol,vector<string> grid, int row, int col, map <int, vector<int>> shapeMapX, map <int, vector<int>> shapeMapY);
tuple<map <int, vector<int>>, map <int, vector<int>>, vector<string>> iterRow (int direction,int node, int corrRow, int corrCol,vector<string> grid, int row, int col, map <int, vector<int>> shapeMapX, map <int, vector<int>> shapeMapY);
pair<map <int, vector<int>>, map <int, vector<int>>> sectionShapes(vector<string> grid, map <int, vector<int>> shapeMapX , map <int, vector<int>> shapeMapY);
void findEdges(float **distances, bool **edges, int **closestCoordsX, int **closestCoordsY, float minDistance, map <int, vector<int>> shapeMapX , map <int, vector<int>> shapeMapY);
pair<vector<float> , vector<float> > findCenters(map <int, vector<int>> shapeMapX , map <int, vector<int>> shapeMapY);
vector<int> findConflict(bool **edges, vector<float> centersX);
pair<map <int, vector<int>>, map <int, vector<int>>> nodeSplitting(int **closestCoordsX, int **closestCoordsY,int node1,int node2,int node3,int vertX,int vertY , map <int, vector<int>> shapeMapX , map <int, vector<int>> shapeMapY);
map <string, vector<int>> findNodesToSplit(vector<int> confNodes, bool **edges, float minDistance, int **closestCoordsX, int **closestCoordsY, map <int, vector<int>> shapeMapX, map <int, vector<int>> shapeMapY);
void plotConflicts(int **closestCoordsX, int**closestCoordsY, int numOfVert, vector<float> centersX, vector<float> centersY, bool **edges);
pair<vector<int>, vector<string>> assignMask(bool **edges, vector<string> grid, map <int, vector<int>> shapeMapX , map <int, vector<int>> shapeMapY);
tuple<map <int, vector<int>>, map <int, vector<int>>, bool>  iterateLoop (float minDistance, map <int, vector<int>> shapeMapX , map <int, vector<int>> shapeMapY);

int main()
{
    // extract layout from file and store in grid
    string line;
    vector<string> grid;
    ifstream layoutFile;
    layoutFile.open("layout.txt");

    if (!layoutFile) 
    {
        cerr << "Unable to open file datafile.txt";
        return 0;   // call system to stop
    }else{
        while (getline(layoutFile,line)) {
            grid.push_back(line);
        }
    }

    layoutFile.close();

    // allow user to change min distance

    string choice;
    float minDistance;

    cout << "\nDefault min distance = 2.80 blocks\n";
    cout << "If you would like to change it, enter y \n";
    cout << "Otherwise, enter anything\n";

    if(getline(cin, choice))
    {
        if(choice =="")
        {
            minDistance = 2.80;
        } else if(choice == "y")
        {
            cin >> minDistance;
        }
        else {
            minDistance = 2.80;
        }
    }


    //////////////////////////////////////////////////////////////////////////
    // GRAPH CONSTRUCTION 
    // SECTION THE SHAPES INTO NODES WITH COORDINATES STORED IN VECTOR MAPS

    map <int, vector<int>> shapeMapX;
    map <int, vector<int>> shapeMapY;

    // function that identifies touching blocks and groups them as nodes
    tie(shapeMapX,shapeMapY) = sectionShapes(grid,shapeMapX,shapeMapY);

    // CALCULATE SHORTEST DISTANCES TO FIND EDGES 

    // define some variables
    int numOfVert = shapeMapX.size();   // number of nodes
    float **distances;                  // matrix of smallest distance between each node
    distances = new float *[numOfVert];
    bool **edges;                       // matrix identifying edge conflicts between nodes (0 = no conflict, 1 = conflict)
    edges = new bool *[numOfVert];
    int **closestCoordsX;               // holds X coordinate of closest distance point between each node
    closestCoordsX = new int *[numOfVert];
    int **closestCoordsY;               // holds Y coordinate of closest distance point between each node
    closestCoordsY = new int *[numOfVert];
    vector<float> centersX;             // holds Y coordinate of shape center point
    vector<float> centersY;             // holds X coordinate of shape center point
    vector<int> confNodes;              // holds nodes identified as in a conflict loop
    // initialize array values
    for (int row = 0; row < numOfVert; row++)    
    {
        distances[row] =  new float [numOfVert];
        edges[row] = new bool [numOfVert];
        closestCoordsX[row] = new int [numOfVert];
        closestCoordsY[row] = new int [numOfVert];
        for (int col = 0; col < numOfVert; col++)
        {
            distances[row][col] = numOfVert * 2;
            edges[row][col] = false;
            closestCoordsX[row][col] = 0;
            closestCoordsY[row][col] = 0;
        }
    }

    // measure the smalles distance between blocks and identifies which one have edges that are below the minimum distance
    findEdges(distances, edges, closestCoordsX, closestCoordsY, minDistance, shapeMapX, shapeMapY);

    // CONSTRUCT A GRAPH TO VISUALIZE THE CONFLICTS

    // finds the center points of each block to be used in the graph construction in plotConflicts
    tie(centersX,centersY) = findCenters(shapeMapX, shapeMapY);
    // plots the conflict graph using centers from findCenters
    plotConflicts(closestCoordsX, closestCoordsY, numOfVert,centersX, centersY,edges);

    //////////////////////////////////////////////////////////////////////////
    // CONFLICT CYCLE DETECTION

    // conflict cycle detector, which finds loops with an odd number of vertices
    confNodes = findConflict(edges, centersX);
    // figures out which nodes should be split to resolve conflicts using the shoftest distance between blocks
    map <string, vector<int>> conflicts = findNodesToSplit(confNodes, edges, minDistance, closestCoordsX, closestCoordsY, shapeMapX, shapeMapY);

    // initialize variables to hold new values after first node splitting
    map <int, vector<int>> shapeMapX2;      // holds new X coordinates of feature shapes after first node splitting
    map <int, vector<int>> shapeMapY2;      // holds new Y coordinates of feature shapes after first node splitting
    tuple<map <int, vector<int>>,map <int, vector<int>>,bool> iterationReturn  (shapeMapX, shapeMapY,false);


    //////////////////////////////////////////////////////////////////////////
    // NODE SPLITTING

    for(int c = 0; c < conflicts["node1"].size();c++) // for each solution to conflict
    {
        // splits nodes at location identified by findNodesToSplit
        tie(shapeMapX2,shapeMapY2) = nodeSplitting(closestCoordsX, closestCoordsY,conflicts["node1"][c],conflicts["node2"][c],conflicts["node3"][c],conflicts["vertX"][c],conflicts["vertY"][c],shapeMapX,shapeMapY);
        // iterative loop used when splitting nodes and checking for new conflicts
        iterationReturn = iterateLoop(minDistance,shapeMapX2, shapeMapY2);
        // if this split is a sucess, store it and break out of the cycle
        if(get<2>(iterationReturn) == true)
        {
            shapeMapX2 = get<0>(iterationReturn);
            shapeMapY2 = get<1>(iterationReturn);
            break;
        }
    }

    // if successful, recalculate all the values for the new graph

    if(get<2>(iterationReturn) == true)
    {
        // initialize variables to hold new values after first node splitting
        int numOfVert2 = shapeMapX2.size(); // number of nodes after first splitting
        float **distances2;                 // matrix of smallest distance between each node after splitting
        distances2 = new float *[numOfVert2];
        bool **edges2;                      // matrix identifying edge conflicts between nodes after splitting
        edges2 = new bool *[numOfVert2];
        int **closestCoordsX2;              // holds X coordinate of closest distance point after splitting
        closestCoordsX2 = new int *[numOfVert2];
        int **closestCoordsY2;              // holds Y coordinate of closest distance point after splitting
        closestCoordsY2 = new int *[numOfVert2];
        vector<float> centersX2;            // holds X coordinate of shape center point after splitting
        vector<float> centersY2;            // holds Y coordinate of shape center point after splitting
        vector<int> masks;                  // holds mask assignments
        // initialize array values
        for (int row = 0; row < numOfVert2; row++)    
        {
            distances2[row] =  new float [numOfVert2];
            edges2[row] = new bool [numOfVert2];
            closestCoordsX2[row] = new int [numOfVert2];
            closestCoordsY2[row] = new int [numOfVert2];
            for (int col = 0; col < numOfVert2; col++)
            {
                distances2[row][col] = numOfVert2 * 2;
                edges2[row][col] = false;
                closestCoordsX2[row][col] = 0;
                closestCoordsY2[row][col] = 0;
            }
        }

        // measure the smalles distance between blocks and identifies which one have edges that are below the minimum distance
        findEdges(distances2, edges2, closestCoordsX2, closestCoordsY2, minDistance, shapeMapX2, shapeMapY2);
        // finds the center points of each block to be used in the graph construction in plotConflicts
        tie(centersX2,centersY2) = findCenters(shapeMapX2, shapeMapY2);

        cout << "Pass? " << "YES" << "\n\n";
        // if all conflicts are resolved, assigns each block a mask
        tie(masks,grid)= assignMask(edges2,grid,shapeMapX2,shapeMapY2);
        // plots the conflict graph using centers from findCenters
        plotConflicts(closestCoordsX2, closestCoordsY2, numOfVert2,centersX2, centersY2,edges2);
    } else {
        cout << "Pass? " << "FAIL - cannot clear conflicts \n\n";
    }
    
    ofstream outputFile;
    outputFile.open ("output.txt");

    // print grid
    for (int i = 0; i < grid.size(); i++)
    {
        std::cout << grid[i] << "\n";
        outputFile << grid[i] << "\n";
    }
    cout<<"\n\n";
    outputFile.close();



    // print some values that may be interesting:
/*
    for (int row = 0; row < numOfVert2; row++)    //for each row
    {
        for (int col = 0 ; col < numOfVert2; col++)   //for each column
        {
            cout << edges2[row][col] << " , ";
        }
        cout << '\n';
    }
    cout << "\n\n";

    int shape = shapeMapX.size();
    vector<int> X;
    vector<int> Y;
    
    for(int i = 0 ; i < shape ; i++)
    {
        X = shapeMapX[i];
        Y = shapeMapY[i];

        for (int j = 0 ; j < X.size() ; j++)
        {
            cout << '(' << X[j] << ", " << Y[j] << ") \n";
        }
        cout << "\n";
    }

    cout << "\n\n";

    
    for(int i = 0 ; i < shapeMapX2.size() ; i++)
    {
        for (int j = 0 ; j < shapeMapX2[i].size() ; j++)
        {
            cout << '(' << shapeMapX2[i][j] << ", " << shapeMapY2[i][j] << ") \n";
        }
        cout << "\n";
    }

    cout << "\n\n";

    for (int row = 0; row < shape; row++)    //for each row
    {
        for (int col = 0 ; col < shape; col++)   //for each column
        {
            cout << edges[row][col] << " , ";
        }
        cout << '\n';
    }

    
    for(int i = 0; i < conflicts["node1"].size(); i++)
    {
        cout << '(' << conflicts["node1"][i] << ", ";
        cout << conflicts["node2"][i] << ", ";
        cout << conflicts["node3"][i] << ") \n";
        cout << '(' << conflicts["vertX"][i] << ", ";
        cout << conflicts["vertY"][i] << ") \n\n";
    }

    
    cout << "\n\n";

    for (int i = 0 ; i < centersX.size() ; i++)
    {
        cout << centersX[i] << '\n';
    }

    cout << "\n\n";

    for (int i = 0 ; i < centersY.size() ; i++)
    {
        cout << centersY[i] << '\n';
    }

    for (int i = 0 ; i < centersX.size() ; i++)
    {
        cout << centersX[i] << '\n';
    }


    for (int row = 0; row < confNodes.size(); row++)    //for each row
    {
            cout << confNodes[row]<< ", ";
        
    }
    cout << "\n\n";
*/


    return 0;
}

// identifies the different blocks in the design and stores their coordinates in vectors
pair<map <int, vector<int>>, map <int, vector<int>>> sectionShapes(vector<string> grid, map <int, vector<int>> shapeMapX , map <int, vector<int>> shapeMapY)
{
    tuple<map <int, vector<int>>,map <int, vector<int>>,vector<string> >shapeItReturn (shapeMapX, shapeMapY,grid);
    int node = 0;
    char nextRow = '_';
    char nextCol = '_';
    char nextCol2 = '_';
    int corrCol = 0;
    int corrRow = 0;
    int corrCol2 = 1;

    map <int, vector<int>>::iterator c;

    for (int row = 0; row < grid.size(); row++)    //for each row
    {
        for (int col = 0 ; col < grid[row].size(); col++)   //for each column
        {         
            if (grid[row][col] == 'X')  // check if this section has a shape that has not been labeled yet
            {
                corrCol = 0;
                corrRow = 0;

                nextCol = grid[row][col];
                nextRow = grid[row][col];

                while ((col + corrCol) < grid[row].size() && nextCol == 'X')    // follow the whole column to keep same letter if connected
                {
                    int direction = 1;
                    shapeItReturn = iterRow(direction,node, corrRow, corrCol,grid,row,col,shapeMapX,shapeMapY);
                    shapeMapX = get<0>(shapeItReturn);
                    shapeMapY = get<1>(shapeItReturn);
                    grid = get<2>(shapeItReturn);

                    direction = -1;
                    shapeItReturn = iterRow(direction,node, corrRow, corrCol,grid,row,col,shapeMapX,shapeMapY);
                    shapeMapX = get<0>(shapeItReturn);
                    shapeMapY = get<1>(shapeItReturn);
                    grid = get<2>(shapeItReturn);

                    corrCol++;
                    corrRow = 0;
                    if ((col + corrCol) < grid[row].size())
                    {
                        nextCol = grid[row+corrRow][col+corrCol];
                    }
                    nextRow = nextCol;
                }
                node++;  // next label
            }
        }
    }

    return make_pair(shapeMapX, shapeMapY);
}

// is used to iterate through columns recursively when identifying blocks in sectionShapes
tuple<map <int, vector<int>>, map <int, vector<int>>, vector<string>> iterCol (int node, int corrRow, int corrCol,vector<string> grid, int row, int col, map <int, vector<int>> shapeMapX, map <int, vector<int>> shapeMapY)
{
    tuple<map <int, vector<int>>,map <int, vector<int>>,vector<string>> shapeItReturn (shapeMapX, shapeMapY,grid);
    map <int, vector<int>>::iterator c;
    char nextCol2;
    int corrCol2 = corrCol + 1;
    int corrRow2;

    if((col + corrCol2) < grid[row].size())
    {
        nextCol2 = grid[row+corrRow][col+corrCol2];
    }
    while(nextCol2 == 'X' && (col + corrCol2) < grid[row].size())   // check for connections to the right
    {
        grid[row+corrRow][col+corrCol2] = 'a';

        c = shapeMapX.find(node);
        if (c == shapeMapX.end()) // if key has not been created, add an empty vector
        {
            shapeMapX.insert(pair<int,vector<int> >(node, vector<int>()));
            shapeMapY.insert(pair<int,vector<int> >(node, vector<int>()));
        }
        shapeMapX[node].push_back(col+corrCol2);
        shapeMapY[node].push_back(row+corrRow);

        if (corrCol2 > 0)    // if we are 2 away from original row, start branching off
        {
            int direction = 1;
            corrRow2 = corrRow + 1;
            shapeItReturn = iterRow(direction, node, corrRow2, corrCol2,grid,row,col,shapeMapX,shapeMapY);
            shapeMapX = get<0>(shapeItReturn);
            shapeMapY = get<1>(shapeItReturn);
            grid = get<2>(shapeItReturn);

            direction = -1;
            corrRow2 = corrRow - 1;
            shapeItReturn = iterRow(direction, node, corrRow2, corrCol2,grid,row,col,shapeMapX,shapeMapY);
            shapeMapX = get<0>(shapeItReturn);
            shapeMapY = get<1>(shapeItReturn);
            grid = get<2>(shapeItReturn);
        }

        corrCol2++;

        if((col + corrCol2) < grid[row].size())
        {
            nextCol2 = grid[row+corrRow][col+corrCol2];
        }
    }

    corrCol2 = corrCol - 1;
    if((col + corrCol2)  >= 0)
    {
        nextCol2 = grid[row+corrRow][col+corrCol2];
    }
    while(nextCol2 == 'X' && (col + corrCol2) >= 0)     // check for connection to the left
    {
        grid[row+corrRow][col+corrCol2] = 'a';

        c = shapeMapX.find(node);
        if (c == shapeMapX.end()) // if key has not been created, add an empty vector
        {
            shapeMapX.insert(pair<int,vector<int> >(node, vector<int>()));
            shapeMapY.insert(pair<int,vector<int> >(node, vector<int>()));
        }
        shapeMapX[node].push_back(col+corrCol2);
        shapeMapY[node].push_back(row+corrRow);

        if (corrCol2 > 0)    // if we are 2 away from original row, start branching off
        {
            int direction = 1;
            corrRow2 = corrRow + 1;
            shapeItReturn = iterRow(direction, node, corrRow2, corrCol2,grid,row,col,shapeMapX,shapeMapY);
            shapeMapX = get<0>(shapeItReturn);
            shapeMapY = get<1>(shapeItReturn);
            grid = get<2>(shapeItReturn);

            direction = -1;
            corrRow2 = corrRow - 1;
            shapeItReturn = iterRow(direction, node, corrRow2, corrCol2,grid,row,col,shapeMapX,shapeMapY);
            shapeMapX = get<0>(shapeItReturn);
            shapeMapY = get<1>(shapeItReturn);
            grid = get<2>(shapeItReturn);
        }

        corrCol2--;

        if((col + corrCol2)  >= 0)
        {
            nextCol2 = grid[row+corrRow][col+corrCol2];
        }
    }
    get<0>(shapeItReturn) = shapeMapX;
    get<1>(shapeItReturn) = shapeMapY;
    get<2>(shapeItReturn) = grid;

    return shapeItReturn;

}

// is used to iterate through rows recursively when identifying blocks in sectionShapes
tuple<map <int, vector<int>>, map <int, vector<int>>, vector<string>> iterRow (int direction, int node, int corrRow, int corrCol,vector<string> grid, int row, int col, map <int, vector<int>> shapeMapX, map <int, vector<int>> shapeMapY)
{
    tuple<map <int, vector<int>>,map <int, vector<int>>,vector<string> >shapeItReturn (shapeMapX, shapeMapY,grid);
    char nextRow;
    int corrRow2 = corrRow;
    map <int, vector<int>>::iterator c;
    
    if ((row + corrRow) < grid.size())
    {
        nextRow = grid[row+corrRow][col+corrCol];
    }

    while ((row + corrRow) < grid.size() && nextRow == 'X') // follow the whole row to keep same letter if connected
    {
        grid[row+corrRow][col+corrCol] = 'a';

        c = shapeMapX.find(node);
        if (c == shapeMapX.end()) // if key has not been created, add an empty vector
        {
            shapeMapX.insert(pair<int,vector<int> >(node, vector<int>()));
            shapeMapY.insert(pair<int,vector<int> >(node, vector<int>()));
        }
        shapeMapX[node].push_back(col+corrCol);
        shapeMapY[node].push_back(row+corrRow);

        if (corrRow > -1)    // if we are 2 away from original row, start branching off
        {
            shapeItReturn = iterCol(node, corrRow, corrCol,grid,row,col,shapeMapX,shapeMapY);
            shapeMapX = get<0>(shapeItReturn);
            shapeMapY = get<1>(shapeItReturn);
            grid = get<2>(shapeItReturn);
        }

        if(direction == 1)
        {
            corrRow++;

        } else {
            corrRow--;
        }

        
        if ((row + corrRow) < grid.size())
        {
            nextRow = grid[row+corrRow][col+corrCol];
        }
    }

    get<0>(shapeItReturn) = shapeMapX;
    get<1>(shapeItReturn) = shapeMapY;
    get<2>(shapeItReturn) = grid;

    return shapeItReturn;

}

// measure the smalles distance between blocks and identifies which one have edges that are below the minimum distance
void findEdges(float **distances, bool **edges, int **closestCoordsX, int **closestCoordsY, float minDistance, map <int, vector<int>> shapeMapX , map <int, vector<int>> shapeMapY)
{
    float distCal = 0.0;
    int mapSize = shapeMapX.size();
    int nextNode = 0;
    vector<int> nodeX;
    vector<int> nodeY;
    vector<int> nextNodeX;
    vector<int> nextNodeY;

    for(int node = 0; node <mapSize; node++) //for every shape/node
    {
        for(int nextNode = 0; nextNode <mapSize; nextNode++) //for every shape/node
        {
            if(node == nextNode)
            {
                distances[node][nextNode] = 0;
                distances[nextNode][node] = 0;

            } else {
                nodeX = shapeMapX[node];
                nodeY = shapeMapY[node];
                nextNodeX = shapeMapX[nextNode];
                nextNodeY = shapeMapY[nextNode];

                for(int itPt = 0; itPt < nodeX.size();itPt++) //for every point in first shape
                {
                    for(int itPt2 = 0; itPt2 < nextNodeX.size();itPt2++) //for every point in second shape
                    {
                        distCal = sqrt( pow( abs(nodeY[itPt] - nextNodeY[itPt2]), 2) + pow(abs(nodeX[itPt] - nextNodeX[itPt2]), 2) ) ;

                        if (distances[node][nextNode] > distCal)
                        {
                            distances[node][nextNode] = distCal;
                            distances[nextNode][node] = distCal;
                            if ( distCal < minDistance)
                            {
                                edges[node][nextNode] = true;
                                edges[nextNode][node] = true;
                                closestCoordsX[node][nextNode] = nodeX[itPt];
                                closestCoordsY[node][nextNode] = nodeY[itPt];
                                closestCoordsX[nextNode][node] = nextNodeX[itPt2];
                                closestCoordsY[nextNode][node] = nextNodeY[itPt2];
                            }
                            if(0 == distCal )
                            {
                                edges[node][nextNode] = false;
                                edges[nextNode][node] = false;
                            }

                        }
                    }
                }
            }
        }
    }
}

// conflict cycle detector, which finds loops with an odd number of vertices
vector<int> findConflict(bool **edges, vector<float> centersX)
{
    vector<float> d;
    vector<bool> visited;
    vector<int> confNodes;

    for (int i = 0 ; i < centersX.size() ; i++)
    {
        d.push_back(-numeric_limits<float>::infinity());
        visited.push_back(false);
    }

    queue<int> Q;
    int k = 0;
    int j = 0;

    Q.push(0);
    d[0] = 0;

    bool found = false;
    
    while (!Q.empty())
    {
        j = Q.front();
        Q.pop();

        for (int adjNode = 0; adjNode < centersX.size() ; adjNode++)
        {
            if(edges[j][adjNode] == true)
            {
                k = adjNode;

                if(d[k] >= 0)
                {
                    if(d[k] == d[j])
                    {
                        found = true;
                        break;
                    }
                }
                else
                {
                    d[k] = d[j] + 1;
                    Q.push(k);
                }
            }
        }
        if(found == true)
        {
            break;
        }
    }



    list<int> L;

    map <int, int> F;
    F.insert(pair<int,int>(j, -1));
    F.insert(pair<int,int>(k, -1));

    queue<int> Q2;
    Q2.push(j);
    Q2.push(k);

    int r = 0;
    int s = 0;
    map <int,int>::iterator check;

    if(found == true)
    {
        while (!Q2.empty())
        {
            if(!L.empty())
            {
                break;
            }
            r = Q2.front();
            Q2.pop();

        for (int adjNode = 0; adjNode < centersX.size() ; adjNode++)
            {
                if(edges[r][adjNode] == true)
                {
                    s = adjNode;

                    if((d[s] + 1) == d[r] )
                    {
                        if(visited[s] == true)
                        {
                            int f;
                            L.push_front(s);
                            L.push_back(r);

                            f = F[r];
                            while (f != -1)
                            {
                                L.push_back(f);
                                f = F[f];
                            }
                            f = F[s];
                            while (f != -1)
                            {
                                L.push_back(f);
                                f = F[f];
                            }
                            break;
                        }

                        check = F.find(s);
                        if (check == F.end())
                        {
                            F.insert(pair<int,int>(s, -1));
                        }

                        F[s] = (r);
                        visited[s] = true;
                        Q2.push(s);
                    }

                }
            }
        }

        int fd = 0;

        for(auto i = L.begin();i != L.end(); ++i)
        {
            for(int j = 0; j< confNodes.size(); ++j)
            {
                if(confNodes[j] == *i)
                {
                    fd = 1;
                }
            }
            if(fd == 0)
            {
                confNodes.push_back(*i);
            }
        }

    }

    return confNodes;

}

// figures out which nodes should be split to resolve conflicts using the shoftest distance between blocks
map <string, vector<int>> findNodesToSplit(vector<int> confNodes, bool **edges, float minDistance, int **closestCoordsX, int **closestCoordsY, map <int, vector<int>> shapeMapX, map <int, vector<int>> shapeMapY)
{
    float dist;
    float vertX;
    float vertY;
    int node1;
    int node2;
    int node3;
    int mapInd = 0;

    map <string, vector<int>> conflicts;

    conflicts.insert(pair<string,vector<int> >("node1", vector<int>()));
    conflicts.insert(pair<string,vector<int> >("node2", vector<int>()));
    conflicts.insert(pair<string,vector<int> >("node3", vector<int>()));
    conflicts.insert(pair<string,vector<int> >("vertX", vector<int>()));
    conflicts.insert(pair<string,vector<int> >("vertY", vector<int>()));

    for (int i = 0 ; i < confNodes.size() ; i++)
    {
        for (int j1 = 0 ; j1 < confNodes.size() ; j1++)
        {
            int j2 = j1+1;
            for(; j2 < confNodes.size() ; j2++)
            {
                node1 = confNodes[i];
                node2 = confNodes[j1];
                node3 = confNodes[j2];
                if(edges[node1][node2] ==true && edges[node1][node3] ==true && node2!=node3)
                {

                    int X1 = closestCoordsX[node1][node2];
                    int X2 = closestCoordsX[node1][node3];
                    int Y1 = closestCoordsY[node1][node2];
                    int Y2 = closestCoordsY[node1][node3];

                    vertX = min(X1,X2) + (abs(X1 - X2)/2);
                    vertY = min(Y1,Y2) + (abs(Y1 - Y2)/2);

                    int currentNodeChar = node1;
                    vector<int> xVect = shapeMapX[currentNodeChar];
                    vector<int> yVect = shapeMapY[currentNodeChar];

                    int closX = vertX;
                    int closY = vertY;
                    float distXY = 20;

                    for (int Ind = 0 ; Ind < xVect.size() ; Ind++)
                    {
                        float d = sqrt( pow( abs(yVect[Ind] - vertY), 2) + pow(abs(xVect[Ind] - vertX), 2) );

                        if(X1 != xVect[Ind] || Y1 != yVect[Ind])
                        {  
                            if(X2 != xVect[Ind] || Y2 != yVect[Ind])
                            {    
                                if(d < distXY)
                                {
                                    distXY = d;
                                    closX = xVect[Ind];
                                    closY = yVect[Ind];
                                }
                            }
                        }
                    }

                    if (closX != vertX)
                    {
                        vertX = closX;
                    }
                    if (closY != vertY)
                    {
                        vertY = closY;
                    }

                    int X12 = closestCoordsX[node2][node1];
                    int X22 = closestCoordsX[node3][node1];
                    int Y12 = closestCoordsY[node2][node1];
                    int Y22 = closestCoordsY[node3][node1];

                    float dist1 = sqrt( pow( abs(Y12 - vertY), 2) + pow(abs(X12 - vertX), 2) );
                    float dist2 = sqrt( pow( abs(Y22 - vertY), 2) + pow(abs(X22 - vertX), 2) );
                    dist = min(dist1,dist2) ;

                    if(minDistance<dist)
                    {
                        conflicts["node1"].push_back(node1);
                        conflicts["node2"].push_back(node2);
                        conflicts["node3"].push_back(node3);
                        conflicts["vertX"].push_back(vertX);
                        conflicts["vertY"].push_back(vertY);
                    }
                }

            }
        }
    }
    return conflicts;
}

// splits nodes at location identified by findNodesToSplit
pair<map <int, vector<int>>, map <int, vector<int>>> nodeSplitting(int **closestCoordsX, int **closestCoordsY, int node1,int node2,int node3,int vertX,int vertY , map <int, vector<int>> shapeMapX , map <int, vector<int>> shapeMapY)
{
    cout << "\n************* \n";
    cout << "TRYING: \n";
    cout << "current:" << node1 << "\n";
    cout << "first:" << node2 << "\n";
    cout << "second:" << node3 << "\n";
    cout << "(" << vertX << ", ";
    cout << vertY << ") \n";
    cout << "************* \n\n";

    bool found = false;
    int newNode = shapeMapX.size() ;

    vector<int> originalX  = shapeMapX[node1];
    vector<int> originalY  = shapeMapY[node1];
    
    map <int, vector<int>> eShapeMapX = shapeMapX;
    map <int, vector<int>> eShapeMapY = shapeMapY;

    int closestX1 = closestCoordsX[node1][node2];
    int closestX2 = closestCoordsX[node1][node3];
    int closestY1 = closestCoordsY[node1][node2];
    int closestY2 = closestCoordsY[node1][node3];

    eShapeMapX.insert(pair<int,vector<int> >(newNode, vector<int>()));
    eShapeMapY.insert(pair<int,vector<int> >(newNode, vector<int>()));


    for (int Ind = 0 ; Ind < originalX.size() ; Ind++)
    {

       float d2 = sqrt( pow( abs(originalY[Ind] - closestY1), 2) + pow(abs(originalX[Ind] - closestX1), 2) );
       float d3 = sqrt( pow( abs(originalY[Ind] - closestY2), 2) + pow(abs(originalX[Ind] - closestX2), 2) );

       if(d2 > d3)
       {
            eShapeMapX[newNode].push_back(originalX[Ind]);
            eShapeMapY[newNode].push_back(originalY[Ind]);
            for (int s = 0; eShapeMapX[node1].size() > s; s++)
            {
                if(eShapeMapX[node1][s] == originalX[Ind] && eShapeMapY[node1][s] == originalY[Ind])
                {
                    eShapeMapX[node1].erase(eShapeMapX[node1].begin()+s);
                    eShapeMapY[node1].erase(eShapeMapY[node1].begin()+s);
                    break;

                }

            }

       }

    }

    return make_pair(eShapeMapX, eShapeMapY);
}

// plots the conflict graph using centers from findCenters
void plotConflicts(int **closestCoordsX, int**closestCoordsY, int numOfVert, vector<float> centersX, vector<float> centersY, bool **edges)
{
    vector<int> x1;
    vector<int> y1;
    vector<int> x2;
    vector<int> y2;
    vector<float> xs;
    vector<float> ys;
    //vector<int>::iterator line;
    int line = 0;
    int col = 0;

    /*
    for (int row = 0; row < numOfVert ; ++row)
    {
        col = row + 1;
        for ( ; col < numOfVert ; ++col)
        {
            if(closestCoordsX[row][col] > 0)
            {
                x1.push_back(closestCoordsX[row][col]);
                y1.push_back(closestCoordsY[row][col]);
                x2.push_back(closestCoordsX[col][row]);
                y2.push_back(closestCoordsY[col][row]);
                cout << '(' << x1.back() << ", " << y1.back() << "), \n";
                cout << '(' << x2.back() << ", " << y2.back() << "), \n\n";

            }
        }
    }
    */

    plt::figure();

    for (int row = 0; row < numOfVert ; ++row)
    {
        col = row + 1;
        for ( ; col < numOfVert ; ++col)
        {
            if(edges[row][col] == 1)
            {
                xs = {centersX[row], centersX[col]};
                ys = {-centersY[row], -centersY[col]};
            }
            // xs = {x1[line], x2[line]};
            //ys = {-y1[line], -y2[line]};
            plt::plot(xs,ys);
        }
    }
    
    plt::show();
}

// finds the center points of each block to be used in the graph construction in plotConflicts
pair<vector<float> , vector<float> > findCenters(map <int, vector<int>> shapeMapX , map <int, vector<int>> shapeMapY)
{
    vector<float> centersX;
    vector<float> centersY;
    vector<int> vectX = shapeMapX[0];
    vector<int> vectY = shapeMapY[0];

    for(int i2 = 0 ; i2 < shapeMapX.size(); i2++) // for each shape
    {
        // average each coordinate
        vectX = shapeMapX[i2];
        vectY = shapeMapY[i2];
        centersX.push_back(0);
        centersY.push_back(0);

        for (int i = 0; i < vectX.size() ; ++i)
        {
            centersX[i2] += vectX[i];
            centersY[i2] += vectY[i];
        }

        centersX[i2] /= vectX.size();
        centersY[i2] /= vectY.size();
    }

    return make_pair(centersX , centersY);
}

// if all conflicts are resolved, assigns each block a mask
pair<vector<int>, vector<string>> assignMask(bool **edges, vector<string> grid, map <int, vector<int>> shapeMapX , map <int, vector<int>> shapeMapY)
{
    vector<int> masks(shapeMapX.size(),3);
    int color = 1;
    queue<int> Q;
    int shape = 0;
    int x = 0;
    int y = 0;
    int i = 0;

    Q.push(0);
    
    while (!Q.empty())
    {
        shape = Q.front();
        Q.pop();

        if (shape == 0){
            masks[0] = 0;
        }

        for(int e = 0; e < masks.size(); ++e)
        {
            if(edges[shape][e] == 1 && masks[e]==3)
            {

                if(masks[shape] == 0)
                {
                    masks[e] = 1;
                } else {
                    masks[e] = 0;
                }

                Q.push(e);
            }
        }
    }

    for(int i = 0; i < shapeMapX.size(); ++i)
    {
        for(int j = 0; j < shapeMapX[i].size(); ++j)
        {
            x = shapeMapX[i][j];
            y = shapeMapY[i][j];

            char c;

            if(masks[i] == 3)
            {
                 c = '0' + 0;
            }else{
                 c = '0' + masks[i];
            }

            grid[y][x] = c;
        }
    }

    return make_pair(masks,grid);
}

// iterative loop used when splitting nodes and checking for new conflicts
tuple<map <int, vector<int>>, map <int, vector<int>>, bool>  iterateLoop (float minDistance, map <int, vector<int>> shapeMapX , map <int, vector<int>> shapeMapY)
{
    int numOfVert = shapeMapX.size();
    float **distances;         // smallest distance between each shape
    distances = new float *[numOfVert];
    bool **edges;
    edges = new bool *[numOfVert];
    int **closestCoordsX;    // holds X coordinate of closest distance point
    closestCoordsX = new int *[numOfVert];
    int **closestCoordsY;    // holds Y coordinate of closest distance point
    closestCoordsY = new int *[numOfVert];

    for (int row = 0; row < numOfVert; row++)    
    {
        distances[row] =  new float [numOfVert];
        edges[row] = new bool [numOfVert];
        closestCoordsX[row] = new int [numOfVert];
        closestCoordsY[row] = new int [numOfVert];
        for (int col = 0; col < numOfVert; col++)
        {
            distances[row][col] = numOfVert * 2;
            edges[row][col] = false;
            closestCoordsX[row][col] = 0;
            closestCoordsY[row][col] = 0;
        }
    }


    findEdges(distances, edges, closestCoordsX, closestCoordsY, minDistance, shapeMapX, shapeMapY);

    vector<float> centersX;
    vector<float> centersY;

    tie(centersX, centersY) = findCenters(shapeMapX, shapeMapY);

    //plotConflicts(closestCoordsX, closestCoordsY, numOfVert,centersX, centersY,edges);

    tuple<map <int, vector<int>>,map <int, vector<int>>,bool> iterationReturn (shapeMapX, shapeMapY,false);

    vector<int> confNodes = findConflict(edges, centersX);
    if(confNodes.size() == 0)
    {
        get<2>(iterationReturn) = true;
        return iterationReturn; 
    }

    map <string, vector<int>> conflicts = findNodesToSplit(confNodes, edges, minDistance, closestCoordsX, closestCoordsY, shapeMapX, shapeMapY);
    if(conflicts["node1"].size() == 0)
    {
        return iterationReturn;
    }

    for(int c = 0; c < conflicts["node1"].size();c++) // for each solution to conflict
    {
        tie(shapeMapX,shapeMapY) = nodeSplitting(closestCoordsX, closestCoordsY,conflicts["node1"][c],conflicts["node2"][c],conflicts["node3"][c],conflicts["vertX"][c],conflicts["vertY"][c],shapeMapX,shapeMapY);

        iterationReturn = iterateLoop(minDistance,shapeMapX, shapeMapY);
        
        if(get<2>(iterationReturn) == true)
        {
            return iterationReturn;
        }
    }
    return iterationReturn; 
    
}

#include <iostream>
#include <math.h>
#include "matplotlibcpp.h"
#include <vector>
using namespace std;
namespace plt = matplotlibcpp;

char sectionShapes(string *grid, int gridSize)
{
    char let = 'a';
    char nextRow = '_';
    char nextCol = '_';
    char nextCol2 = '_';
    int corrCol = 0;
    int corrRow = 0;
    int corrCol2 = 1;

    for (int row = 0; row < gridSize; row++)    //for each row
    {
        for (int col = 0 ; col < gridSize; col++)   //for each column
        {         

            if (grid[row][col] == 'X')  // check if this section has a shape that has not been labeled yet
            {
                corrCol = 0;
                corrRow = 0;

                nextCol = grid[row][col];
                nextRow = grid[row][col];

                while ((col + corrCol) < gridSize && nextCol == 'X')    // follow the whole column to keep same letter if connected
                {
                    while ((row + corrRow) < gridSize && nextRow == 'X') // follow the whole row to keep same letter if connected
                    {
                        grid[row+corrRow][col+corrCol] = let;   // replace the X with the current label
                        corrRow++;

                        if ((row + corrRow) < gridSize)
                        {
                            nextRow = grid[row+corrRow][col+corrCol];
                        }

                        if (corrRow > 1)    // if we are 2 away from original row, start branching off
                        {
                            corrCol2 = corrCol + 1;
                            if((col + corrCol2) < gridSize)
                            {
                                nextCol2 = grid[row+corrRow][col+corrCol2];
                            }
                            while(nextCol2 == 'X' && (col + corrCol2) < gridSize)   // check for connections to the right
                            {
                                grid[row+corrRow][col+corrCol2] = let;
                                corrCol2++;

                                if((col + corrCol2) < gridSize)
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
                                grid[row+corrRow][col+corrCol2] = let;
                                corrCol2++;

                                if((col + corrCol2)  >= 0)
                                {
                                    nextCol2 = grid[row+corrRow][col+corrCol2];
                                }
                            }
                        }
                    }

                    corrCol++;
                    corrRow = 0;
                    if ((col + corrCol) < gridSize)
                    {
                        nextCol = grid[row+corrRow][col+corrCol];
                    }
                    nextRow = nextCol;
                }
                let++;  // next label
            }
        }
    }

    return let;
}

pair<map <char, vector<int>>, map <char, vector<int>>>findEdges(string *grid, float **distances, bool **edges, int **closestCoordsX, int **closestCoordsY, int gridSize, float minDistance, map <char, vector<int>> shapeMapX , map <char, vector<int>> shapeMapY)
{
    int corrCol = 0;
    int corrRow = 0;
    float distCal = 0.0;
    char currentChar = '_';
    char otherChar = '_';
    char a = 'a';
    int currentInd = 0;
    int otherCharInd = 0;
    map <char, vector<int>>::iterator c;

    for (int row = 0; row < gridSize; row++)    //for each row
    {
        for (int col = 0 ; col < gridSize; col++)   //for each column
        {
            currentChar = grid[row][col];

            if (currentChar != '_')
            {
                c = shapeMapX.find(currentChar);
                if (c == shapeMapX.end()) // if key has not been created, add an empty vector
                {
                    shapeMapX.insert(pair<int,vector<int> >(currentChar, vector<int>()));
                    shapeMapY.insert(pair<int,vector<int> >(currentChar, vector<int>()));
                }
                shapeMapX[currentChar].push_back(col);
                shapeMapY[currentChar].push_back(row);
                currentInd = currentChar - a;

                for ( corrRow = 0; corrRow < gridSize ; corrRow++ )
                {
                    for ( corrCol = 0; corrCol < gridSize ; corrCol++ )
                    {
                        otherChar = grid[corrRow][corrCol];
                        otherCharInd = otherChar - a;

                        if (otherChar != '_')
                        {
                            if (currentChar != otherChar)
                            {
                                distCal = sqrt( pow( abs(row - corrRow), 2) + pow(abs(col - corrCol), 2) ) ;

                                if (distances[currentInd][otherCharInd] > distCal)
                                {
                                    distances[currentInd][otherCharInd] = distCal;
                                    if (distCal < minDistance)
                                    {
                                        edges[currentInd][otherCharInd] = true;
                                        closestCoordsX[currentInd][otherCharInd] = col;
                                        closestCoordsY[currentInd][otherCharInd] = row;
                                    }

                                }
                            } else {
                                distances[currentInd][otherCharInd] = 0;
                            }
                            
                        }

                    }

                }


            }

        }
    }

    return make_pair(shapeMapX, shapeMapY);

}

pair<vector<float> , vector<float> > findCenters(map <char, vector<int>> shapeMapX , map <char, vector<int>> shapeMapY, vector<float> centersX, vector<float> centersY)
{
    char k = 'a';
    char end  = shapeMapX.count(k);
    int i = 0;
    int i2 = 0;

    vector<int> vectX;
    vector<int> vectY;

    for( ; end==1 ; ) // for each shape
    {
        // average each coordinate
        i2 = k - 'a';
        vectX = shapeMapX[k];
        vectY = shapeMapY[k];
        centersX.push_back(0);
        centersY.push_back(0);

        for (i = 0; i < vectX.size() ; ++i)
        {
            centersX[i2] += vectX[i];
            centersY[i2] += vectY[i];
        }

        centersX[i2] /= vectX.size();
        centersY[i2] /= vectY.size();

        k++;
        end  = shapeMapX.count(k);
    }

    return make_pair(centersX , centersY);

}

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


int main()
{
    int gridSize = 13;
    float minDistance = 2.80;

    string grid[] = {
        "__XXX_XXX____",
        "_____________",
        "_XXXX_XXX_XXX",
        "____________X",
        "XXXX________X",
        "X___________X",
        "X_XXXXXXXX_XX",
        "__X__________",
        "__X____X_____",
        "__X____X_____",
        "_____________",
        "__XXXXX______",
        "_____________"
    };

    char let = 'a';

    ////////////////////////////////////////////////////////////////
    // GRAPH CONSTRUCTION

    // SECTION THE SHAPES USING LETTERS TO IDENTIFY THEM
    let = sectionShapes(grid,gridSize);

    ////////////////////////////////////////////////////////////////////////////
    // CALCULATE SHORTEST DISTANCES TO FIND EDGES 

    int numOfVert = let - 'a';
    int vertCoords [numOfVert][2];                  // midpoint of each shape
    float **distances;         // smallest distance between each shape
    distances = new float *[numOfVert];
    bool **edges;
    edges = new bool *[numOfVert];
    int **closestCoordsX;    // holds X coordinate of closest distance point
    closestCoordsX = new int *[numOfVert];
    int **closestCoordsY;    // holds Y coordinate of closest distance point
    closestCoordsY = new int *[numOfVert];
    map <char, vector<int>> shapeMapX;
    map <char, vector<int>> shapeMapY;
    vector<float> centersX; // holds Y coordinate of shape center point
    vector<float> centersY; // holds X coordinate of shape center point


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

    tie(shapeMapX,shapeMapY) = findEdges(grid, distances, edges, closestCoordsX, closestCoordsY, gridSize, minDistance, shapeMapX, shapeMapY);
    tie(centersX,centersY) = findCenters(shapeMapX, shapeMapY, centersX, centersY);

    plotConflicts(closestCoordsX, closestCoordsY, numOfVert,centersX, centersY,edges);



    for (int i = 0; i < (sizeof(grid)/sizeof(grid[0])); i++)
    {
        std::cout << grid[i] << "\n";
    }


    for (int row = 0; row < numOfVert; row++)    //for each row
    {
        for (int col = 0 ; col < numOfVert; col++)   //for each column
        {
            cout << edges[row][col] << " , ";
        }
        cout << '\n';
    }

    for (int i = 0 ; i < centersX.size() ; i++)
    {
        cout << centersX[i] << '\n';
    }
    

    return 0;
}

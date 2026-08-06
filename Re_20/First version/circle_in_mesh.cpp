//the last code was...
//also a huge sucess but basu sir said that it would be better if i dont use that F, S, I differenciation becasue...
//..i am just wasting memory ! but rater save interface nodes in one array, the solid interior nodes in one array...
//...this should saves us much more memory and time!  

//the last code was...
//..a huge sucess BUT i need to find the percentage of the cell that is inside the solid and the percentage of the cell..
//that is outside! and rather than representing it as S,F, I i think i should represent it as 1 ,0, and something...
//...in between for the interface! 

//!!!!!----CIRCLE-------!!!!!

//here i am gonna try and find where the cells completley of fluid, completley of solid and the interface are located! 
//the main idea for the circle is going to be the distance formula! the idea is gonna be the mesh is gonna be a class and...
//..what i will do is gonna initiall the class mesh will be set as F = 1, S = 0, I = 0, the idea is that 1 is yes and 0 is..
//..No. the idea is gonna be initaly everything will be set as fluid and then we will check if the node is inside, outside..
//..or on the interface somewhere...

//also what do i know? 
//- center of the circle
//- radius of the circle

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <array>
#include <fstream>
#include <chrono>
#include <sstream> 

#include <iomanip> //this if for only printing a particular number of digits for exaple printing 0.4 when the actual number is 0.44332 or somthing

using namespace std;


class property
{
    public:
        int x;
        int y;
        double ratio;

        property()
        {
            ratio = 0;
        }
        

};


int main()
{

    auto star = chrono::high_resolution_clock::now();

    int nodes_x = 51; //both of these are gonna be the x and y length of the reagion we will be searching!
    int nodes_y = 51;

    int total_nodes = nodes_x * nodes_y;

    int center_x = 26; //center location of the circle!
    int center_y = 26;

    //hmm the biggest problem i think i will face is that my code is completly written in a 1D array! all the while the center..
    //...of the circle is given as some x and y, also for properly using the distance formula i need the x and y! 
    //..maybe a fucntion to turn the 1D array to 2D? 
    // I have already done this! something i like int x = i % nodes_x, and int y = i / nodes_x woudl work i guess!

    double radius = 10;

    vector<property> mesh(total_nodes);

    vector<property> interface;
    vector<property> solid;

    for(int i = 0; i < total_nodes; i++)
    {
        int x = i % nodes_x;
        int y = i / nodes_x;
        mesh[i].x = x;
        mesh[i].y = y;

    }

    for(int i = 0; i < total_nodes; i++)
    {
        if( mesh[i].y < center_y - (radius + 1) || mesh[i].x < center_x - (radius + 1) || mesh[i].y > center_y + (radius + 1) || mesh[i].x > center_x + (radius + 1))
        {
            continue;
        }
        else if(mesh[i].y > center_y - (radius + 1) && mesh[i].x > center_x - (radius + 1) && mesh[i].y < center_y + (radius + 1) && mesh[i].x < center_x + (radius + 1))
        {
            double dist = sqrt( pow((mesh[i].x - center_x), 2) + pow((mesh[i].y - center_y), 2));

            if(dist > radius + 0.5)
            {
                continue;
            }
            else if(dist < radius - 0.5)
            {
                solid.push_back(mesh[i]);
                mesh[i].ratio = 1.0;
            }

            else if(dist <= radius + 0.5 && dist >= radius - 0.5)
            {
                interface.push_back(mesh[i]);
            }
        }
    }

    //i have already found the interfaces BUT i want to know what percentage of the cell is filled by solid!
    //i already have an idea on how to do it! the cell must be divided again into more parts, for now 5 x 5 parts..
    //..but how will i do it?
    //the first idea that came into my mind is why dont i just creat more nodes there like at 12, 15 location...
    //i can add stuff like 11.5, 11.75, 12, 12.25, 12.5 now those are 5 pices! 

    

    
        
        // double solid_count = 0.0;
        // double fluid_count = 0.0;

        // int div = 5;

        // double val = 1.0 / div;

        // int x = i % nodes_x;
        // int y = i / nodes_x;

        // double x_new = x - 0.5;
        // double y_new = y - 0.5;
        // double x_high = x + 0.5;
        // double y_high = y + 0.5;


        for(int j = 0; j < interface.size(); j++)
        {
            double solid_count = 0.0;
            double fluid_count = 0.0;

            int div = 5;

            double val = 1.0 / div;

            int x = interface[j].x;
            int y = interface[j].y;

            int mesh_id = (y * nodes_x) + x;

            double x_new = x - 0.5;
            double y_new = y - 0.5;
            double x_high = x + 0.5;
            double y_high = y + 0.5;

            for(int a = 0; a < div ; a++)
            {
                for(int b = 0; b < div; b++)
                {
                    double x_current = x_new + (a*val) + (val/2);
                    double y_current = y_new + (b * val) + (val/2);
                    double dist = sqrt( pow((x_current - center_x), 2) + pow((y_current - center_y), 2));

                    if(dist <= radius)
                    {
                        solid_count++;
                    }
                    
                    
                }
            }
            mesh[mesh_id].ratio = solid_count/(div*div);

        }

    


    for(int i = 0; i < total_nodes; i++)
    {
        if(i % nodes_x == 0)
        {
            cout << "\n";
        }
        
        cout << fixed << setprecision(1) << mesh[i].ratio;
        cout << "";

    }

    ofstream outDom("Domain.csv");

    outDom << "x,y,Domain\n";


    
    for(int i = 0; i < mesh.size(); i++) 
    {

        outDom << mesh[i].x << "," << mesh[i].y << "," << mesh[i].ratio << "\n";

    }

  

    outDom.close();

    cout << "Domain saved as Domain.csv!" << endl;

    //i also want to plot some of the locations on cricle 

    
    double x_start = center_x - radius;
    
    //the number of nodes needed in circle is set using this cir_res

    int cir_res = 20;

    vector<double> circle_x(cir_res + 1);
    vector<double> circle_y(cir_res + 1);

    double theta_div = (2.0 * M_PI) / cir_res;

    for(int i = 0; i < (cir_res); i++)
    {
        double theta = i * theta_div;
        
        circle_x[i] = center_x + (radius * cos(theta));
        circle_y[i] = center_y + (radius * sin(theta));

    }

    ofstream outCir("Circle.csv");

    outCir << "x,y\n";


    
    for(int i = 0; i < cir_res ; i++) 
    {

        outCir << circle_x[i] << "," << circle_y[i] << "\n";

    }

    outCir << circle_x[0] << "," << circle_y[0] << "\n";

    outCir.close();

    cout << "Domain saved as Circle.csv!" << endl;


}
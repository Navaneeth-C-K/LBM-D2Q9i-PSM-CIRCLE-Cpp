//------------------------The beginning on Paritially saturated method------------------------------------
//ok let's try and underysand the idea behind PSM and why it is used! 
//lets go for why first! 
//The problem is rn i will be implimenting the PSM on the square cylinder! this will be..
//..a test case so that i can unserstadn if the solution i am getting in both are correct or not!
//BUT what is the actual use? 
//the actual use comes when we are looking at non grid allogned geometries! for example a circle!
//when putting a cirlce in a grid the problem is that we wont be able to have all the nodes on the..
//..circle (obviously!) so we use PSM! the idea is to set fluids as 0 and solids as 1..
//..and the interface cells as something between 0 and 1! 
//Cells! this the word cell will sound intresting as we have been dealing with nodes so far!..
//..the idea is 
//...         *  ___ *    
//              | . |
//              |___|   this is a (terrble, i would say horrible) representation of a cell!
//            *       *
//a cell has a node at it's center and is surrounderd on all side by other nodes a cell here is..
//..also 1 lattice unit in size.
//An intresting question here would be why is there a need for a cell! why not just use nodes? 
//well this is where the 0 , 1 and between 0 and 1 comes into play! 
//a cell with no solid inside will be 0,  a cell with only solid inside will be 1 and a cell...
//..with let's say 60% solid and 40% fluid will be 0.6! this is the idea behind the PSM
//then we use an equation which will have some kind of factor which will help the solver understand..
//..what is solid, fluid and interface! easy peasy!

//------------------------The beginning of the code---------------------------------------

//The code is:
//-> D2Q9i
//-> BGK
//-> error calculation using L2 norm
//-> Drag and lift calculation using momentum exhchange method
//-> PSM (implimeting)
//-> previous output as input!


//problem statement: Flow past a square cylinder using PSM
//d is gonna be the size of the square cylinder in lattice units! every other distace is gonna be multiples of it

//Outputs:
//Velocity at each node, density at each node, drag vs time and lift vs time all saved as csv files!
//other outputs that are printed out are the Breer et al comparison value


#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <array>
#include <fstream>
#include <chrono>
#include <sstream> //this is for reading the output files!

using namespace std;

class Data
{
    public:
        double x;
        double y;
        double psi;

};

class vector2D
{
    
    public:
        double x, y;
        vector2D()
        {
            x = 0;
            y = 0;
        }

        vector2D(double x_, double y_)
        {
            x = x_;
            y = y_;
        }

    double operator*(const vector2D &rhs)
    {
        return (this->x*rhs.x + this->y*rhs.y);
    }

    vector2D operator*(double val)
    {
        return vector2D(this->x*val, this->y*val);
    }

    vector2D operator+(const vector2D &rhs)
    {
        return vector2D(this->x + rhs.x, this->y + rhs.y);
    }

    vector2D operator=(const vector2D &rhs)
    {
        return vector2D(x = rhs.x, y = rhs.y);
    }

    vector2D operator/(double val)
    {
        if(val==0){
        cout << "Not defined" << endl;
        return vector2D(this->x, this->y);
        }
        else{
        return vector2D(this->x/val, this->y/val);
        }
    }
void print_vector()
    {
        cout << "X = " << x << ", Y = " << y <<endl;
    }

    friend ostream& operator<<(ostream &out, const vector2D &v);

    
};

ostream& operator<<(ostream &out, const vector2D &v)
    {
        out << "X = " << v.x << ", Y = " << v.y <<endl;

        return out;
    }

class direction
{
    //array<vector2D, 9> xi;

    public:
    array<vector2D, 9> xi;
        direction()
        {
            xi[0] = vector2D(0, 0);
            xi[1] = vector2D(1, 0);
            xi[2] = vector2D(0, 1);
            xi[3] = vector2D(-1, 0);
            xi[4] = vector2D(0, -1);
            xi[5] = vector2D(1, 1);
            xi[6] = vector2D(-1, 1);
            xi[7] = vector2D(-1, -1);
            xi[8] = vector2D(1, -1);
        }

        vector2D operator[](int rhs) const
        {
            return xi[rhs];
        }

};

class weight
{
    array<double, 9> w;

    public:
        weight()
        {
            w[0] = 4.0/9.0;
            w[1] = 1.0/9.0;
            w[2] = 1.0/9.0;
            w[3] = 1.0/9.0;
            w[4] = 1.0/9.0;
            w[5] = 1.0/36.0;
            w[6] = 1.0/36.0;
            w[7] = 1.0/36.0;
            w[8] = 1.0/36.0;


        }

        weight operator*(const double &rhs) 
        {
            weight product;
            for(int i = 0; i<9; i++)
            {
                product.w[i] =  w[i]*rhs;
            }
            return product;
        }

        double operator[](int rhs) const
        {
            return w[rhs];
        }
};

class node_property

{
    
    weight W;
    //direction Xi; will make this private later!
    
 

    public:
        direction Xi;
        array<double, 9> f;
        double rho;
        vector2D vel;
        double ratio;
        int x;
        int y;

        

        void equilibrium(double rho, vector2D u, double Cs, double rho_0)//here i want to multiply w[i] with a variable rho 
        {

            for(int i= 0; i<9;i++)
            {
                f[i] =  (W[i])*(rho + ( rho_0 * ((u*Xi[i])/(Cs*Cs)) + (0.5* (u*Xi[i])*(u*Xi[i])/(Cs*Cs*Cs*Cs)) - (0.5*(u*u)/(Cs*Cs)) ));
            }

        }

        double equilibrium(double rho, vector2D u, double Cs, int n, double rho_0) //this will help me create direction specific f_equilibrium (as it is impo for boundary node calculations)
        {
            //node_property equi; 
                // ratio = 0.0;
                return  (W[n])*(rho + ( rho_0 * ((u*Xi[n])/(Cs*Cs)) + (0.5* (u*Xi[n])*(u*Xi[n])/(Cs*Cs*Cs*Cs)) - (0.5*(u*u)/(Cs*Cs)) ));
            
        }



};


int main()
{

//All the setup stuff happens here!
    auto start = chrono::high_resolution_clock::now();

    //all the set up stuff happens here! 

    int d = 32; 

    int radius = d/2;

    int center_x = 10*d; 
    int center_y = 8*d;

    double Re = 40;
    int nodes_y = 16*d;
    int nodes_x = 30*d; 
    int pad_x = nodes_x + 2;//this i think will be much more convenient way  of adding those graveyard nodes.. padding!!
    int pad_y = nodes_y + 2;
    int total_nodes = pad_x * pad_y;//here i am just finding the total nodes


    //this is for finding wake length!
    int cylinder_end_x = center_x + radius;

    int centerline_y = center_y;

    int wake_start_idx = (centerline_y * pad_x) + cylinder_end_x + 1;


    double L_x = 20;// unit meter
    double L_y = 1;//unit meter 



    double Ma = 0.1;
    //double Cs = 1/sqrt(3.0);
    double Cs = 0.5773502691896257;


    double U_star = Ma * Cs;
    double neu_star = (U_star * d) / Re; //why nodes_y?? easy remmeber equation for reynolds number! nodes_y is just the h distance between the two plates!
    //double neu_star = (1.5 * U_star * d) / Re;// i made a huge stupid mistake here! forgot to multiply it U_star with 1.5 as that is U_max!

    double Tau = (3*neu_star) + 0.5;


    double density = 1000;
    double rho = 1;

    double rho_0 = 1;

    double rho_constant = rho;
    double rho_in = rho;


    double u_in = 1;
    vector2D vel_zero(0, 0);


    double total_time = 600.0;


    double node_length = L_x / nodes_x;
    double time_step = (U_star * node_length) / u_in;
    

    cout << time_step << endl;
    cout << U_star << endl;

    vector<node_property> mesh(total_nodes);

    vector<node_property> temporary(total_nodes);

    vector<node_property> equilibrium_array(total_nodes);

    double drag_final1;
    double drag_final2;
    double drag_final3;
    double drag_final4;

    double lift_final1;
    double lift_final2;
    double lift_final3;
    double lift_final4;






    ifstream inVel("velocity_field_sq_cy.csv");

    if (inVel.is_open())
    {
        cout << "Previous velocity file found!" << endl;

        //Step 1: i am just gonna initialize 
        for(int i = 0; i < total_nodes; i++)
        {
            mesh[i].rho = rho; // rho is 1.0
            mesh[i].vel = vel_zero;
            mesh[i].ratio = 0.0;
        }

        string line;
        getline(inVel, line); //this is just used to ignore the top header lines i.e the x, y, u, v written on top!

        while (getline(inVel, line)) 
        {
            stringstream ss(line);
            string val;
            
            int x, y;
            double u_scaled, v_scaled;

            // remember to always skip commas ck!
            getline(ss, val, ','); x = stoi(val);
            getline(ss, val, ','); y = stoi(val);
            getline(ss, val, ','); u_scaled = stod(val);
            getline(ss, val, ','); v_scaled = stod(val);

            // calculaing the 1D index
            int i = (y * pad_x) + x;

            // when i saved the u and v into my velocity feild earlier i divided both the u and v with u star here i am gonna...
            //..multiply it back to get the accurate value.
            mesh[i].vel.x = u_scaled * U_star;
            mesh[i].vel.y = v_scaled * U_star;
        }

        inVel.close();




        
    }

    ifstream inDens("Density_sq_cy.csv");

    if (inDens.is_open())
    {
        cout << "Previous Density file found!" << endl;

        

        string line;
        getline(inDens, line); //this is just used to ignore the top header lines i.e the x, y, u, v written on top!

        while (getline(inDens, line)) 
        {
            stringstream ss(line);
            string val;
            
            int x, y;
            double dens;

            // remember to always skip commas ck
            getline(ss, val, ','); x = stoi(val);
            getline(ss, val, ','); y = stoi(val);
            getline(ss, val, ','); dens = stod(val);
            

            // calculaing the 1D index
            int i = (y * pad_x) + x;

            // when i saved the u and v into my velocity feild earlier i divided both the u and v with u star here i am gonna...
            //..multiply it back to get the accurate value.
            mesh[i].rho = dens ;
            
        }

        inDens.close();

        //setting up the equilibrium again.
        for(int i = 0; i < total_nodes; i++) 
        {
            mesh[i].equilibrium(mesh[i].rho, mesh[i].vel, Cs, rho_0);
        }

    }

    // else
    // {
    //     for(int i = 0; i < total_nodes; i++)
    //     {
    //         mesh[i].rho = rho;
    //         mesh[i].vel = vel_zero;
    //         mesh[i].ratio = 0.0;
    //         mesh[i].equilibrium(mesh[i].rho, mesh[i].vel, Cs, rho_0);
    //     }
    // }


    else
    {
        // Moved the variable declaration INSIDE the else block
        vector2D vel_init(U_star, 0.0); 

        for(int i = 0; i < total_nodes; i++)
        {
            mesh[i].rho = rho;
            mesh[i].vel = vel_init; // Initialize to U_star
            mesh[i].ratio = 0.0;
            
            // This initializes the populations f_i to equilibrium based on U_star
            mesh[i].equilibrium(mesh[i].rho, mesh[i].vel, Cs, rho_0); 
        }
    }

    

    //******************************************************************************************************************** */

    int t = 0;
    double total_time_lu =  total_time/time_step;

    //*******************************************************************************************

    //the time steps start from here! 

    vector<double>  u_array(total_nodes);
    
    double saved_squared_sum = 0;

    double error = 1;

    int steps = 0;

    int count_out1;
    int count_out2;

    vector<double> lift_out1;
    vector<double> lift_out2;

    vector<double> drag_saved;
    vector<double> lift_saved;

//----------circle stuff----------------------------------------
    
    vector<node_property> interface;
    vector<node_property> solid;

    vector2D vel_solid(0,0);

    for(int i = 0; i < total_nodes; i++)
    {
        int x = i % pad_x;
        int y = i / pad_x;
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

    for(int j = 0; j < interface.size(); j++)
    {
        double solid_count = 0.0;
        double fluid_count = 0.0;

        int div = 5;

        double val = 1.0 / div;

        int x = interface[j].x;
        int y = interface[j].y;

        int mesh_id = (y * pad_x) + x;

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

    double C_d_last_time = 0.0;

//---------cricle stuff complete------------------------------------


    ofstream outDrag("Drag_Circle.dat");

    ofstream outLift("Lift_circle.dat");

    ofstream outRecir("Recirculation_circle.dat");

//------------while loop starts HERE!!---------------------------------------

    //there is a huge disclamer that is needed here! the problem is that i have tried for like 800,000..
    //..time stpes and it still did not converge using the L2 scheme for velocity!
    //so i am gonna keep the while loop as it is and then when Cd gets a acuracy of something like 0.0001..
    //..then we will break out of the while loop!
    //while(error > 0.001) //this error will be for the Cd term!!
    while(t< total_time_lu)
    {

        double squared_sum = 0.0;
        double sum_square_diff = 0.0;

        double Force = 0.0;

        double drag1 = 0.0; 
        double drag2 = 0.0; 
        double drag3 = 0.0; 
        double drag4 = 0.0; 

        double lift1 = 0.0;
        double lift2 = 0.0;
        double lift3 = 0.0;
        double lift4 = 0.0;

        vector<double> lifts1; //debugging of lift purpose
        vector<double> lifts2; //debugging of lift purpose
        vector<int> loc1; // debugging of lift purpose
        vector<int> loc2; // debugging of lift purpose

        int count1 = 0; //debugging of lift purpose
        int count2 = 0; //debugging of lift purpose
        
        

        //the colission step!

       for(int i = nodes_x+3; i <= total_nodes-(nodes_x+4); i++) //colission step
        {
            for(int j = 0; j < 9 ; j++)
            {
                if(i % (nodes_x + 2) == 0 || (i + 1) % (nodes_x + 2) == 0) {continue;}

                else
                {
                    double Equi_solid = mesh[i].equilibrium(mesh[i].rho, vel_solid , Cs, j, rho_0);

                    double N_d =  mesh[i].equilibrium(mesh[i].rho , mesh[i].vel, Cs, j, rho_0);
                    //temporary[i].f[j] = mesh[i].f[j] - ((mesh[i].f[j]-N_d)/Tau); // here i am doing the collision and saving it as temporary array
                    //the main difference in PSM will come here! 
                    //here i will be using the PSM equation to calculate the temporatry f[j] at i
                    //but the ratio of solid to fluid in the interface will be needed before hand! So i need to..
                    //..find it before hand and save it as mesh[i].ratio
                    //also i think is it fine to have the cell and node have the same name as each cell only...
                    //..contains one node! BUT the problem i have with it is that the boundary..
                    //..specifically at the top and bottom walls of the fluid flow domain is made of nodes! 
                    //..but again i dont think i have to worry about that for now! 
                    //BUT this is something i most likely will have to come back to again!(i think)
                    //BUT then again PSM can actually be used at the top and bottom walls too right?
                    //but why not just use zho he bounce back there and PSM at object in the flow? 
                    //is there a problem in that? something to think about!

                    double Omega_BGK = -(mesh[i].f[j] - N_d)/Tau;

                    double Omega_S = (Equi_solid) - (mesh[i].f[j]) + ((1.0 - (1.0/Tau))*(mesh[i].f[j] - N_d));

                    double B = ((mesh[i].ratio * (Tau - 0.5)) / ((1 - mesh[i].ratio) + (Tau - 0.5)));

                    temporary[i].f[j] = mesh[i].f[j] + ((1.0 - B)*Omega_BGK) + (B* Omega_S);
                    
                    
                    drag1= drag1 - (B* Omega_S * mesh[i].Xi[j].x );
                    lift1= lift1 - (B * Omega_S * mesh[i].Xi[j].y );
                   
                }

            }    
        }
    

        //streaming step!


        for(int i = nodes_x+3; i <= total_nodes-(nodes_x+4); i++) //streaming step
        {
            int x = i % pad_x;
            int y = i / pad_x;

            if(i % (pad_x) == 0 || (i + 1) % (pad_x) == 0) {continue;}

            
            

            else
            {

                for(int j = 0; j < 9 ; j++)
                {
            
                    int a;
                    if(j == 0){a = 0;}
                    if(j == 1){a = 1;}
                    if(j == 2){a = nodes_x+2;}
                    if(j == 3){a = -1;}
                    if(j == 4){a = -nodes_x-2;}
                    if(j == 5){a = nodes_x +2 + 1;}
                    if(j == 6){a = nodes_x + 2 - 1;}
                    if(j == 7){a = -nodes_x - 2 - 1;}
                    if(j == 8){a = -nodes_x - 2 + 1;}

                    int target = i + a;

                    mesh[target].f[j] = temporary[i].f[j];     

                }  
            }
        }


        //************************************************************************** */
        //Boundary setups!!

        
        int top_left = total_nodes - (nodes_x + nodes_x + 3);
        int top_right = total_nodes - (nodes_x + 4);
        int bot_left = nodes_x + 3;
        int bot_right = nodes_x + 2 + nodes_x;

        int sq_top_left;
        int sq_bottom_left;
        int sq_top_right;
        int sq_bottom_right;

        int k = 0; // this variable i just made to set the inlet velocti in a parabolic form just think of this a somekind of count!
            //as we know the we will start from the bottom left corner and move up so within this else if the k must increase like the the distance in y axis (hope that makes sense future navaneeth ck)

        for(int i = nodes_x+3; i <= total_nodes-(nodes_x+4); i++)
        {
            


            if(i % (pad_x) == 0 || (i + 1) % (pad_x) == 0) {continue;}

            

            //TOP WALL
            else if(i<(total_nodes-(nodes_x + 4)) && i >(total_nodes-(nodes_x + nodes_x + 3)))
            {
                // int a = 0;
                // int b = 1;

                // double u_x = 0;
                // double u_y = 0;

                // vector2D wall_vel(u_x, u_y); 
             
                // double f_equal = mesh[i].f[0] + mesh[i].f[1] + mesh[i].f[3];
                // double f_plus = mesh[i].f[6] + mesh[i].f[2] + mesh[i].f[5];
            
                // //mesh[i].rho = (f_equal + (2*f_plus))/(1+ (a*u_x) + (b*u_y));
                // mesh[i].rho = (f_equal + (2*f_plus))-( (a*u_x) + (b*u_y));

                // mesh[i].f[4] = mesh[i].f[2] - mesh[i].equilibrium(mesh[i].rho, wall_vel, Cs, 2, rho_0) + mesh[i].equilibrium(mesh[i].rho, wall_vel, Cs,  4, rho_0);

                // mesh[i].f[7] = ((mesh[i].rho*(-u_x - u_y)) + mesh[i].f[1] + mesh[i].f[2] - mesh[i].f[3] - mesh[i].f[4] + (2*mesh[i].f[5]))*0.5;
                // mesh[i].f[8] = ((mesh[i].rho*(u_x - u_y)) - mesh[i].f[1] + mesh[i].f[2] + mesh[i].f[3] - mesh[i].f[4] + (2*mesh[i].f[6]))*0.5;

                // mesh[i].vel.x = 0.0;
                // mesh[i].vel.y = 0.0;

                for(int jf =0; jf < 9; jf++)
                {
                    mesh[i].f[jf] = mesh[i - pad_x].f[jf];
                    mesh[i].rho = mesh[i - pad_x].rho;
                    mesh[i].vel.x = mesh[i - pad_x].vel.x;
                    mesh[i].vel.y = 0.0;
                }
            

            }


            //BOTTOM WALL
            else if(i>(nodes_x + 3) && i<(nodes_x + 2 + nodes_x))
            {


                // int a = 0;
                // int b = -1;

                // double u_x = 0;
                // double u_y = 0;

                // vector2D wall_vel(u_x, u_y); 

                // double f_equal = mesh[i].f[3] + mesh[i].f[0] + mesh[i].f[1];
                // double f_plus = mesh[i].f[7] + mesh[i].f[4] + mesh[i].f[8];
            
                // //mesh[i].rho = (f_equal + (2*f_plus))/(1+ (a*u_x) + (b*u_y));
                // mesh[i].rho = (f_equal + (2*f_plus))-( (a*u_x) + (b*u_y));

                // mesh[i].f[2] = mesh[i].f[4] - mesh[i].equilibrium(mesh[i].rho, wall_vel, Cs, 4, rho_0) + mesh[i].equilibrium(mesh[i].rho, wall_vel, Cs,  2, rho_0);

                // mesh[i].f[5] = ((mesh[i].rho*(u_x + u_y)) - mesh[i].f[1] - mesh[i].f[2] + mesh[i].f[3] + mesh[i].f[4] + (2*mesh[i].f[7]))*0.5;
                // mesh[i].f[6] = ((mesh[i].rho*(-u_x + u_y)) + mesh[i].f[1] - mesh[i].f[2] - mesh[i].f[3] + mesh[i].f[4] + (2*mesh[i].f[8]))*0.5;

                // mesh[i].vel.x = 0.0;
                // mesh[i].vel.y = 0.0;

                for(int jf =0; jf < 9; jf++)
                {
                    mesh[i].f[jf] = mesh[i + pad_x].f[jf];
                    mesh[i].rho = mesh[i + pad_x].rho;
                    mesh[i].vel.x = mesh[i + pad_x].vel.x;
                    mesh[i].vel.y = 0.0;
                }

            
            }

            
            //LEFT WALL


           //explanation for variable k given on top!
            else if((i-1)%(pad_x) == 0)//left wall
            {

                int a = -1;
                int b = 0;

                double f_equal = mesh[i].f[0] + mesh[i].f[2] + mesh[i].f[4];
                double f_plus = mesh[i].f[6] + mesh[i].f[3] + mesh[i].f[7];

                double u_y = 0; //as in plane pouisslei flow there is no velocity in y direction!

                // // 1. Get the physical Y coordinate of the current node
                // int y = i / pad_x;

                // // 2. The true physical center of the fluid channel (e.g., node 136.5 for d=34)
                // double physical_center_y = (nodes_y + 1) / 2.0; 

                // // 3. The distance from the center to the wall
                // double max_radius = (nodes_y - 1) / 2.0; 

                // // 4. Calculate the symmetric parabola based strictly on Y coordinates
                // //double u_x = (1.5 * U_star) * (1.0 - pow((y - physical_center_y) / max_radius, 2));
                // double u_x = (1.5 * U_star) * (1.0 - ((y - physical_center_y) / max_radius)*((y - physical_center_y) / max_radius));



                double u_x = U_star; 

                //-------------------------------------------------------------------------------------------------------------------------------
                //our major difference will come here! as we cant just use u_x = U_star as our inlet! are going to set the flow as fully developed right at the inlet itslef
                //so our inlet will be of a parabolic velocity profile with the peak as 1.5x U_star
                //so for this we will have to set each of the inlet nodes to be of a particular velocity each some factor of U_star


                // double center_y = nodes_y / 2.0;
                // double u_x = (1.5 * U_star) * (1.0 - pow((k - center_y) / center_y, 2));


            
                vector2D wall_vel(u_x, u_y); 

                //mesh[i].rho = (f_equal + (2*f_plus))/(1+ (a*u_x) + (b*u_y));
                mesh[i].rho = (f_equal + (2*f_plus))-( (a*u_x) + (b*u_y));

                mesh[i].f[1] = mesh[i].f[3] + ((2.0/3.0) * rho_0 * u_x);

                // mesh[i].f[5] = mesh[i].f[7] - (0.5 * (mesh[i].f[2] - mesh[i].f[4])) + ( (1.0/6.0) * mesh[i].rho*u_x);

                // mesh[i].f[8] = mesh[i].f[6] + (0.5 * (mesh[i].f[2] - mesh[i].f[4])) + ( (1.0/6.0) * mesh[i].rho*u_x);

                mesh[i].f[5] = mesh[i].f[7] - (0.5 * (mesh[i].f[2] - mesh[i].f[4])) + ( (1.0/6.0) * rho_0*u_x);

                mesh[i].f[8] = mesh[i].f[6] + (0.5 * (mesh[i].f[2] - mesh[i].f[4])) + ( (1.0/6.0) * rho_0*u_x);

                mesh[i].vel.x = u_x;
                mesh[i].vel.y = u_y;

                k++;


            }

            //RIGHT WALL
            else if((i+2) % pad_x == 0)
            {
                mesh[i].rho = 1.0;
                //double u_y = mesh[i-1].vel.y; //this was an error! this is only first order accurate! 
                double u_y = ((4.0*mesh[i-1].vel.y) - (mesh[i-2].vel.y))/3.0; //second order accurate
                //double u_x = -1.0 + (mesh[i].f[0] + mesh[i].f[2] + mesh[i].f[4] + 2.0 * (mesh[i].f[1] + mesh[i].f[5] + mesh[i].f[8])) / mesh[i].rho;
                double u_x =  (mesh[i].f[0] + mesh[i].f[2] + mesh[i].f[4] + 2.0 * (mesh[i].f[1] + mesh[i].f[5] + mesh[i].f[8])) - mesh[i].rho;
                vector2D wall_vel(u_x, u_y);
                //mesh[i].f[3] = mesh[i].f[1] - (2.0 / 3.0) * mesh[i].rho * u_x;
                mesh[i].f[3] = mesh[i].f[1] - mesh[i].equilibrium(mesh[i].rho, wall_vel, Cs, 1, rho_0) + mesh[i].equilibrium(mesh[i].rho, wall_vel, Cs,  3, rho_0);

                // mesh[i].f[6] = mesh[i].f[8] - 0.5 * (mesh[i].f[2] - mesh[i].f[4]) - (1.0 / 6.0) * mesh[i].rho * u_x + 0.5 * mesh[i].rho * u_y;
                // mesh[i].f[7] = mesh[i].f[5] + 0.5 * (mesh[i].f[2] - mesh[i].f[4]) - (1.0 / 6.0) * mesh[i].rho * u_x - 0.5 * mesh[i].rho * u_y;
                
                // mesh[i].f[6] = ((mesh[i].rho*(u_y - u_x)) + mesh[i].f[1] - mesh[i].f[2] - mesh[i].f[3] + mesh[i].f[4] + (2*mesh[i].f[8]))*0.5;

                // mesh[i].f[7] = ((mesh[i].rho*(-u_x-u_y)) + mesh[i].f[1] + mesh[i].f[2] - mesh[i].f[3] - mesh[i].f[4] + (2*mesh[i].f[5]))*0.5;

                mesh[i].f[6] = ((rho_0*(u_y - u_x)) + mesh[i].f[1] - mesh[i].f[2] - mesh[i].f[3] + mesh[i].f[4] + (2*mesh[i].f[8]))*0.5;

                mesh[i].f[7] = ((rho_0*(-u_x-u_y)) + mesh[i].f[1] + mesh[i].f[2] - mesh[i].f[3] - mesh[i].f[4] + (2*mesh[i].f[5]))*0.5;



                mesh[i].vel.x = u_x; mesh[i].vel.y = u_y;
            }

         
            
        }



        //i am just gonna copy the values of the nodes right next to the node in the x direction as the f values.
        mesh[top_left] = mesh[top_left + 1];   
        mesh[top_right] = mesh[top_right - 1]; 
        mesh[bot_left] = mesh[bot_left + 1];   
        mesh[bot_right] = mesh[bot_right - 1];

        
        // mesh[sq_top_left] = mesh[sq_top_left + 1];
        // mesh[sq_top_right] = mesh[sq_top_right - 1];
        // mesh[sq_bottom_left] = mesh[sq_bottom_left + 1];
        // mesh[sq_bottom_right] = mesh[sq_bottom_right - 1];

        //for implimenting mometum exchage i will need to need opposites of f... i.e opposite of f1 will be f3 ( -> opposite is <-)
        int opposite[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6}; //so this guy we will be using in the drag and lift calculation! 


        for(int i = nodes_x+3; i <= total_nodes-(nodes_x+4); i++)
        {
            if(i % (nodes_x + 2) == 0 || (i + 1) % (nodes_x + 2) == 0) {continue;} 

           
           
            
            vector2D sum_f_xi_temp(0, 0); 
        
            mesh[i].rho = 0;

            for(int j = 0; j < 9 ; j++)//this loop will just be there for velocity calculation!
            {
                mesh[i].rho = mesh[i].rho + mesh[i].f[j]; // i initially set the special loop just for rho calculation here but i guess i could just calculate the rho when streaming!
                
                sum_f_xi_temp = sum_f_xi_temp + (mesh[i].Xi[j] * mesh[i].f[j]);
            }

            mesh[i].vel = sum_f_xi_temp; //here the divided by rho is removed by me as for d2q9 it is not required (well it is divided by 1..)

                //if(steps <= 10000) {continue;}

            if(steps > 1000 && steps % 100 == 0)
            {
                        
                sum_square_diff = sum_square_diff + pow((mesh[i].vel.x - u_array[i]), 2);

                squared_sum = squared_sum + pow((mesh[i].vel.x), 2);
            }

            u_array[i] = mesh[i].vel.x;
               
            

        
        }


        cout << mesh[nodes_x+3].rho << endl;
        cout << mesh[nodes_x+3+nodes_x].rho << endl;
        cout << mesh[total_nodes-nodes_x-3].rho <<endl;
        cout << mesh[total_nodes -nodes_x-nodes_x-2].rho <<endl;
        cout <<mesh[35].rho <<endl;

        double sum=0;

        //i am trying to make the code run faster so i am getting rid of this!
        // for(int i = nodes_x+3; i < total_nodes-(nodes_x+4); i++)
        // {
        //     sum = sum + mesh[i].rho;
        // }

        cout <<"*********"<<endl;
        cout << sum << endl;
        cout<<"runnig!" <<endl;
        
        cout << t << "/" << total_time_lu << endl;


        // I have commented this out due to the reason i explained right  before the starting of the while loop
        // if(steps > 1000 && steps % 100 == 0)
        // {
        //     error = sqrt(sum_square_diff / squared_sum);
        // } 

        t++;
        steps++;

        saved_squared_sum  = squared_sum;

        drag_final1 = drag1;
       

        lift_final1 = lift1;
        
        

        lift_out1 = lifts1;
        lift_out2 = lifts2;

        //drag_saved.push_back((drag1 + drag2 + drag3 + drag4)/(0.5 * (1.5 * U_star) * (1.5 * U_star) * d));
        //lift_saved.push_back((lift1 + lift2 + lift3 + lift4)/(0.5 * (1.5 * U_star) * (1.5 * U_star) * d));

        
        outDrag << t << "\t" << ((drag1)/(0.5 * ( U_star) * ( U_star) * d)) << "\n";

        outLift << t << "\t" << ((lift1)/(0.5 * (U_star) * (U_star) * d)) << "\n";


        error = ((drag1)/(0.5 * ( U_star) * ( U_star) * d)) - C_d_last_time;


        C_d_last_time = ((drag1)/(0.5 * ( U_star) * ( U_star) * d));

    }

//--------------Post while loop things!!--------------------------------------------

    double time_taken = t*time_step;

    cout << "Time taken for convergence =" << time_taken << endl;

    ofstream outField("velocity_field_sq_cy.csv");

    outField << "x,y,u,v\n";

   

   

    for(int y = 1; y <= nodes_y; y++) 
    {

        for(int x = 1; x <= nodes_x; x++) 
        {

            

            int i = (y * pad_x) + x;

            outField << x << "," << y << "," << mesh[i].vel.x/U_star << "," << mesh[i].vel.y/U_star << "\n";
            //outField << x << "," << y << "," << mesh[i].vel.x << "," << mesh[i].vel.y << "density =" << mesh[i].rho << "\n";


        }

    }

    outField.close();

    cout << "2D velocity field saved to velocity_field_sq_cy.csv!" << endl;

    int gad  = (pad_y + 1) / 2;  

    int gad2 = (gad * pad_x) + nodes_x;

    int gad3 = (gad * pad_x) + 1;


    //outputting desity file

    ofstream outDensity("Density_sq_cy.csv");

    outDensity << "x,y,rho\n";


    // Looping only over the actual fluid domain (1 to nodes)
    for(int y = 1; y <= nodes_y; y++) 
    {

        for(int x = 1; x <= nodes_x; x++) 
        {

            

            int i = (y * pad_x) + x;

            outDensity << x << "," << y << "," << mesh[i].rho << "\n";
            //outField << x << "," << y << "," << mesh[i].vel.x << "," << mesh[i].vel.y << "density =" << mesh[i].rho << "\n";


        }

    }

  

    outDensity.close();

    cout << "2D velocity field saved to Density_sq_cy.csv!" << endl;



    //---------------------



    cout << mesh[gad2-5].vel.x/U_star << endl;
    for(int i = 0 ;i < nodes_x; i++)
    {
        cout << mesh[gad2-i].vel.x/U_star << endl;
    }
    cout << "density*************" << endl;
    for(int i = 0 ;i < nodes_x; i++)
    {
        cout << mesh[gad2-i].rho << endl;
    }
    cout << mesh[gad2].vel.x/U_star << "****" <<  endl;

    cout << mesh[gad3].vel.x/U_star <<endl;

    int dist;

    // for(int i = 0 ;i < nodes_x; i++)
    // {
    //     if(mesh[gad2-i].vel.x < 0)
    //     {
    //         dist = i;
    //     }
    // }

    // cout << dist <<endl;

    


    int wake_length_nodes = 0;
    double neg_vel = 0.0;
    double length;
    for(int i = 0; i < (nodes_x - cylinder_end_x); i++)
    {
        int current_node = wake_start_idx + i;

        // If the velocity is negative, the wake is still going
        if(mesh[wake_start_idx + i].vel.x < 0) 
        {
            wake_length_nodes++;
            neg_vel = mesh[current_node].vel.x;
        }
        else // As soon as velocity hits 0 or positive, the wake has closed!
        {
            double pos_vel = mesh[current_node].vel.x;

            length = ((0-neg_vel)/(pos_vel - neg_vel)) + (wake_length_nodes);
            break; 
        }
    }

    //-------------------calculation of drag and lift coefficient--------------------------------------------

    double coeff_drag1 = drag_final1/(0.5 * ( U_star) * ( U_star) * d);
    


    double coeff_lift1 = lift_final1/(0.5 * (U_star) * ( U_star) * d);

    //----------------End of calculation of drag and lift coefficint-----------------------------------------

    //****************************************************************************************************** 

    //----------------Printing into terminal-----------------------------------------------------------------

    cout << "For Reynolds Number of " << Re << "and d = " << d <<  endl;
    cout << "Wake length in nodes: " << wake_length_nodes << endl;
    cout << "Wake ends at X coordinate: " << cylinder_end_x + wake_length_nodes << endl;
    cout << "Non-dimensional recirculation length : " << length / d << endl;

    

    double total_drag_coeff = coeff_drag1 ;
    double total_lift_coeff = coeff_lift1 ;

    cout << "*********************************" << endl;

    cout << "Totoal drag coefficient = " << total_drag_coeff << endl;
    cout << "Total lift coefficient =  " << total_lift_coeff << endl;




//*****************lift and drag file closing**************************** */

    outDrag.close();

    cout << "Drag saved to Drag_circle.csv!" << endl;



    outLift.close();

    cout << "Lift saved to Lift_circle.csv!" << endl;

    outRecir.close();

    cout << "recirculation length saved to Recirculation_circle.csv" << endl;


    //***********************Lift and drag saved************************* */

    auto end = chrono::high_resolution_clock::now(); 

    chrono::duration<double> elapsed = end - start; //calculation of time taken! 

    cout << "Time taken for the full simulation : " << elapsed.count() << "seconds\n";



}




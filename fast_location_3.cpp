#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Random.h>

#include <vector>
#include <cassert>
#include <iostream>
#include <string>
#include <boost/math/constants/constants.hpp>
#include <fstream>
#include <boost/algorithm/string.hpp>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_3<K, CGAL::Fast_location> Delaunay;
typedef Delaunay::Point Point;
typedef CGAL::Vector_3<K> Vector;
typedef std::vector<Point> br;
typedef boost::tuple<double,br*> R_br;
typedef Delaunay::Cell_handle Cell_handle;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Locate_type    Locate_type;
using namespace std;
int main()
{
	string filename="c:/out/cube.txt";
	std::cout<<"start reading file"<<"["<<filename<<"]"<<endl;

	ifstream ifs(filename);

	if(ifs.fail()){
		std::cout<<"File does not exist."<<endl;
		exit(0);
	}
	std::vector<R_br> data;
	br* branch;
	string line;
	int N;
	double R;

	// read file
	while(getline(ifs,line))
	{
		std::vector<string> words;
		boost::algorithm::split( words, line, boost::algorithm::is_any_of(","));
		//cout << boost::algorithm::join( words, "-" ) << endl;
		char _prefix[20];
		sscanf(words[0].data(),"%s",_prefix);
		string prefix=string(_prefix);
		if(prefix=="P")
		{
			branch=new br();
			sscanf(words[1].data(),"%d",&N);
		}
		if(prefix=="R")
		{
			sscanf(words[1].data(),"%lf",&R);
			//cout<<"new branch N:"<<N<<" R:"<<R<<endl;
			data.push_back(boost::make_tuple(R,branch));
		}
		if(prefix=="C"){
			double x,y,z;
			sscanf(words[1].data(),"%lf",&x);
			sscanf(words[2].data(),"%lf",&y);
			sscanf(words[3].data(),"%lf",&z);
			branch->push_back(Point(x,y,z));
		}
	}

	//Generating points along branches.
	std::vector<Point> exterior;
	std::vector<Point> interior;
	int C=0;
	double pi=boost::math::constants::pi<double>();
	for(vector<R_br>::iterator itr=data.begin();itr!=data.end();itr++,C++)
	{
		double Radius;
		br* _branch;
		boost::tie(Radius,_branch)=*itr;
		for(br::iterator anitr=_branch->begin();anitr!=_branch->end()-1;anitr++)
		{
			Point P=*anitr;
			Point Q=*(anitr+1);
			Vector V=Q-P;
			V=V/std::sqrt(V.squared_length());
			Vector Z=Vector(0,0,1);
			Vector W=CGAL::cross_product<K>(Z,V);
			Vector T=CGAL::cross_product<K>(Z,W);
			exterior.push_back(P);
			for(int i=0;i<12;i++)
			{
				double theta=(double)i/12.*2.*pi;
				Vector N=Radius*(W*std::cos(theta)+T*std::sin(theta));
				Point D=Point(N.x()+0.5*(P.x()+Q.x()),N.y()+0.5*(P.y()+Q.y()),N.z()+0.5*(P.z()+Q.z()));
				exterior.push_back(D);
				for(double scale=0.1;scale<1.0;scale+=0.1)
				{
					Vector N=scale*Radius*(W*std::cos(theta)+T*std::sin(theta));
					Point D=Point(N.x()+0.5*(P.x()+Q.x()),N.y()+0.5*(P.y()+Q.y()),N.z()+0.5*(P.z()+Q.z()));
					interior.push_back(D);				
				}
			}
		}
		exterior.push_back(*_branch->end());
	}
	
	//File write
	/*filename="c:/out/cube.out";
	cout<<"start writing file"<<"["<<filename<<"]"<<endl;

	ofstream ofs(filename);
	ofs<<exterior.size()<<endl;
	for(std::vector<Point>::iterator itr=exterior.begin();itr!=exterior.end();itr++)
	{
		Point P=*itr;
		ofs<<P.x()<<" "<<P.y()<<" "<<P.z()<<endl;
	}
	ofs.close();
	*/
	// building their Delaunay triangulation.
	std::cout<<"start triangulation"<<endl;
	//Delaunay T(exterior.begin(), exterior.end());
	Delaunay T;
	N=0;
	int NN=(exterior.size()/20);
	for(std::vector<Point>::iterator itr=exterior.begin();itr!=exterior.end();itr++,N++)
	{
		T.insert(*itr);
		if(((int)N/NN)*NN==N)
		{
			std::cout<<"*";
		}
	}
	std::cout<<endl;
	std::cout<<"end triangulation"<<endl;
	//int* contains=new int[T.number_of_cells()];
	Locate_type lt;
	int li, lj;
	N=0;
	std::map<Delaunay::Cell_iterator,int> index;
	std::cout<<"start cell search"<<endl;
	
	for(Delaunay::Cell_iterator itr=T.cells_begin();itr!=T.cells_end();itr++,N++)
	{
		index.insert(pair<const Delaunay::Cell_iterator,int>(itr,N));
	}
	Cell_handle before=Cell_handle();
	Cell_handle c=Cell_handle();
	N=0;
	NN=(interior.size()/20);
	for(vector<Point>::iterator itr=interior.begin();itr!=interior.end();itr++,N++)
	{
		if(N==0)
			c = T.locate(*itr, lt,li,lj);
		else
			c=T.locate(*itr,lt,li,lj,before);
		std::map<Delaunay::Cell_iterator,int>::iterator it_c=index.find(c);
		int d=it_c->second;
		before=c;
		if(((int)N/NN)*NN==N)
		{
			std::cout<<"*";
		}
	}
	std::cout<<endl;
	std::cout<<"end cell search"<<endl;
	std::cout << "Press Return To Exit...";
	std::cin.get();
	// performing nearest vertex queries to a series of random points,
	// which is a case where the Fast_location policy is beneficial.
/*	for (int i=0; i<10000; ++i)
		T.nearest_vertex(Point(CGAL::default_random.get_double(0, 20),
		CGAL::default_random.get_double(0, 20),
		CGAL::default_random.get_double(0, 20)));
	*/
	//release memory
	for(vector<R_br>::iterator itr=data.begin();itr!=data.end();itr++)
	{
		R_br a=*itr;
		const br* _branch=boost::get<1>(a);
		delete(_branch);
	}



  return 0;
}

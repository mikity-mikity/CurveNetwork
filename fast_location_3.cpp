#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/Random.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <vector>
#include <cassert>
#include <iostream>
#include <string>
#include <boost/math/constants/constants.hpp>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include<algorithm>	
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned,K>    Vb;
typedef CGAL::Triangulation_data_structure_3<Vb>                    Tds;
//typedef CGAL::Delaunay_triangulation_3<K, Tds,CGAL::Fast_location> Delaunay;
typedef CGAL::Surface_mesh_default_triangulation_3 Delaunay;
typedef	CGAL::Complex_2_in_triangulation_3<Delaunay> C2t3;

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
	string NAME="test";
	string filename="c:/out/"+NAME+".txt";
	std::cout<<"start reading file"<<"["<<filename<<"]"<<endl;

	ifstream ifs(filename);

	if(ifs.fail()){
		std::cout<<"File does not exist."<<endl;
		exit(0);
	}
	std::vector<R_br> data;
	br* branch;
	string line;
	unsigned N;
	double R;

	// read file
	while(getline(ifs,line))
	{
		std::vector<string> words;
		boost::algorithm::split( words, line, boost::algorithm::is_any_of(","));
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
	int num=0;
	double pi=boost::math::constants::pi<double>();
	for(vector<R_br>::iterator itr=data.begin();itr!=data.end();itr++,C++)
	{
		double Radius;
		br* _branch;
		boost::tie(Radius,_branch)=*itr;
		exterior.push_back(*(_branch->begin()));
		for(br::iterator anitr=_branch->begin();anitr!=_branch->end()-1;anitr++)
		{
			Point P=*anitr;
			Point Q=*(anitr+1);
			Vector V=Q-P;
			V=V/std::sqrt(V.squared_length());
			Vector Z=Vector(0,0,1);
			double t=V*Z;
			if(std::fabs(t)>0.9)
			{
				Z=Vector(1,0,0);			
			}
			Vector W=CGAL::cross_product<K>(Z,V);
			Vector T=CGAL::cross_product<K>(V,W);
			//exterior.push_back(P);
			for(int i=0;i<12;i++)
			{
				double theta=(double)i/12.*2.*pi;
				Vector N=Radius*(W*std::cos(theta)+T*std::sin(theta));
				Point D=Point(N.x()+0.5*(P.x()+Q.x()),N.y()+0.5*(P.y()+Q.y()),N.z()+0.5*(P.z()+Q.z()));
				exterior.push_back(D);
			}
		}
		exterior.push_back(*(_branch->end()-1));
	}
	
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
			double t=V*Z;
			if(std::fabs(t)>0.9)
			{
				Z=Vector(1,0,0);			
			}
			Vector W=CGAL::cross_product<K>(Z,V);
			Vector T=CGAL::cross_product<K>(V,W);
			//exterior.push_back(P);
			for(int i=0;i<12;i++)
			{
				double theta=(double)((double)i+0.5)/12.*2.*pi;
				for(double scale=0.5;scale<1.0;scale+=0.05)
				{
					Vector N=scale*Radius*(W*std::cos(theta)+T*std::sin(theta));
					Point D=Point(N.x()+(P.x()),N.y()+(P.y()),N.z()+(P.z()));
					interior.push_back(D);				
				}
			}
		}
	}

	// building their Delaunay triangulation.
	std::cout<<"start triangulation"<<endl;
	
	Delaunay T;
	C2t3 ct(T);
	std::cout<<"exterior.size()"<<exterior.size()<<endl;
	N=0;
	int NN=(exterior.size()/20);
	for(int i=0;i<20;i++)
	{
		std::vector<Point>::iterator be=exterior.begin()+exterior.size()*i/20;
		std::vector<Point>::iterator en=exterior.begin()+exterior.size()*(i+1)/20;
		T.insert(be,en);
		std::cout<<"*";
	}
	std::cout<<endl;
	std::cout<<"end triangulation"<<endl;
	std::cout<<"T.number_of_vertices"<<T.number_of_vertices()<<endl;
	std::cout<<"number_of_cells:"<<T.number_of_cells()<<endl;
	std::cout<<"number_of_finite_cells:"<<T.number_of_finite_cells()<<endl;
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
	std::vector<int> cells;

	for(int i=0;i<T.number_of_cells();i++)
	{
		cells.push_back(0);
	}
	for(vector<Point>::iterator itr=interior.begin();itr!=interior.end();itr++,N++)
	{
		if(N==0)
			c = T.locate(*itr, lt,li,lj);
		else
			c=T.locate(*itr,lt,li,lj,before);
		std::map<Delaunay::Cell_iterator,int>::iterator it_c=index.find(c);
		int d=it_c->second;
		cells[d]++;
		before=c;
		if(((int)N/NN)*NN==N)
		{
			std::cout<<"*";
		}
	}
	std::cout<<endl;
	std::cout<<"end cell search"<<endl;
	cout<<"T.number_of_cells:"<<T.number_of_cells()<<endl;
	cout<<"T.number_of_facets:"<<T.number_of_facets()<<endl;
	cout<<"ct.number_of_facets:"<<ct.number_of_facets()<<endl;

	int max=*(std::max_element(cells.begin(),cells.end()));
	for(int i=0;i<=max;i++)
	{
		std::cout<<i<<":"<<std::count(cells.begin(),cells.end(),i)<<endl;
	}
	
	NN=T.number_of_facets()/20;
	N=0;
	for(Delaunay::Cell_iterator itr=T.cells_begin();itr!=T.cells_end();itr++,N++)
	{
		std::map<Delaunay::Cell_iterator,int>::iterator it_c=index.find(itr);
		int d=it_c->second;
		if(cells[d]>=1)
		{
			ct.add_to_complex(itr,0);
			ct.add_to_complex(itr,1);
			ct.add_to_complex(itr,2);
			ct.add_to_complex(itr,3);
		}
		if(((int)N/NN)*NN==N)cout<<"*";
	}
	cout<<endl;
	cout<<"complex created"<<endl;
	cout<<"T.number_of_facets:"<<T.number_of_facets()<<endl;
	cout<<"ct.number_of_facets:"<<ct.number_of_facets()<<endl;
	std::map<C2t3::Cell_handle,int> space;
	num=0;
	for(C2t3::Facet_iterator itr=ct.facets_begin();itr!=ct.facets_end();itr++,num++)
	{
		space.insert(std::make_pair(itr->first,num));
	}



	//File write

	
	std::map<C2t3::Vertex_handle,int> vIndex;
	filename="c:/out/"+NAME+".out";
	string filename2="c:/out/"+NAME+"_S.out";
	cout<<"start writing file"<<"["<<filename<<"]"<<endl;
	ofstream ofs(filename);
	ofstream ofs2(filename2);
	N=0;
	for(C2t3::Facet_iterator itr=ct.facets_begin();itr!=ct.facets_end();itr++)
	{
		C2t3::Facet opf=ct.opposite_facet(*itr);
		//if(space.find(opf.first)!=space.end())
		{
			for(int i=0;i<4;i++)
			{
				if(i!=itr->second)
				{
					C2t3::Vertex_handle handle=itr->first->vertex(i);
				
					std::map<C2t3::Vertex_handle,int>::iterator pair=vIndex.find(handle);
					if(pair==vIndex.end())
					{
						ofs<<handle->point().x()<<" , "<<handle->point().y()<<" , "<<handle->point().z()<<endl;
						vIndex.insert(std::make_pair(handle,N));
						ofs2<<N<<" ";
						N++;
					}else
					{
						ofs2<<pair->second<<" ";
					}
				}
			}
			ofs2<<endl;
		}
	}
	ofs.close();
	ofs2.close();
	

	std::cout << "Press Return To Exit...";
	std::cin.get();
	// performing nearest vertex queries to a series of random points,
	// which is a case where the Fast_location policy is beneficial.
	//release memory
	for(vector<R_br>::iterator itr=data.begin();itr!=data.end();itr++)
	{
		R_br a=*itr;
		const br* _branch=boost::get<1>(a);
		delete(_branch);
	}
	//delete(cells);



  return 0;
}

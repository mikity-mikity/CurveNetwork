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
typedef CGAL::Delaunay_triangulation_3<K, CGAL::Fast_location> Delaunay;
typedef Delaunay::Point Point;
typedef CGAL::Vector_3<K> Vector;
typedef std::vector<Point> br;
typedef std::vector<Vector> planes;
typedef boost::tuple<double,br*> R_br;
typedef Delaunay::Cell_handle Cell_handle;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Locate_type    Locate_type;
using namespace std;
int main()
{
	string NAME="heart";
	string filename="c:/out/"+NAME+".txt";
	std::cout<<"start reading file"<<"["<<filename<<"]"<<endl;

	ifstream ifs(filename);

	if(ifs.fail()){
		std::cout<<"File does not exist."<<endl;
		exit(0);
	}
	std::vector<R_br> data;
	std::vector<planes*> planeTree;
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
		planes *_planes=new planes();
		planeTree.push_back(_planes);
		for(br::iterator anitr=_branch->begin();anitr!=_branch->end();anitr++)
		{
			Vector V;
			if(anitr==_branch->begin())
			{
				Point P=*anitr;
				Point Q=*(anitr+1);
				V=Q-P;
			}else if(anitr==_branch->end()-1)
			{
				Point P=*(anitr-1);
				Point Q=*(anitr);
				V=Q-P;
			}else
			{
				Point P=*(anitr-1);
				Point Q=*(anitr);
				Point R=*(anitr+1);
				V=(R-P)*0.5;
			}
			V=V/std::sqrt(V.squared_length());
			_planes->push_back(V);
		}
	}
	vector<R_br>::iterator itrA=data.begin();
	vector<planes*>::iterator itrB=planeTree.begin();
	std::cout<<"construct exterior"<<endl;
	while(itrA!=data.end())
	{
		double Radius;
		br* _branch;
		planes* _planes=*itrB;
		boost::tie(Radius,_branch)=*itrA;
		exterior.push_back(*(_branch->begin()));
		br::iterator itrC=_branch->begin();
		planes::iterator itrD=_planes->begin();
		while(itrC!=_branch->end())
		{
			for(int i=0;i<12;i++)
			{
				double theta=(double)i/12.*2.*pi;
				Vector V=*itrD;
				Vector Z=Vector(0,0,1);
				double t=V*Z;
				if(std::fabs(t)>0.9)
				{
					Z=Vector(1,0,0);			
				}
				Vector W=CGAL::cross_product<K>(Z,V);
				Vector T=CGAL::cross_product<K>(V,W);
				W=W/std::sqrt(W.squared_length());
				T=T/std::sqrt(T.squared_length());
				Vector N=Radius*(W*std::cos(theta)+T*std::sin(theta));
				Point D=Point(N.x()+(*itrC).x(),N.y()+(*itrC).y(),N.z()+(*itrC).z());
				exterior.push_back(D);
			}
			itrC++;
			itrD++;
		}
		exterior.push_back(*(_branch->end()-1));		
		itrA++;
		itrB++;
	}
	itrA=data.begin();
	itrB=planeTree.begin();
	std::cout<<"construct interior"<<endl;
	while(itrA!=data.end())
	{
		double Radius;
		br* _branch;
		planes* _planes=*itrB;
		boost::tie(Radius,_branch)=*itrA;
		exterior.push_back(*(_branch->begin()));
		br::iterator itrC=_branch->begin();
		planes::iterator itrD=_planes->begin();
		while(itrC!=_branch->end()-1)
		{
			Vector V1=*itrD;
			Vector V2=*(itrD+1);
			Point P1=*itrC;
			Point P2=*(itrC+1);
			for(double s=0.1;s<1.0;s+=0.2)
			{
				Vector V=V1*s+V2*(1.-s);
				Vector Z=Vector(0,0,1);
				double t=V*Z;
				if(std::fabs(t)>0.9)
				{
					Z=Vector(1,0,0);			
				}
				Vector W=CGAL::cross_product<K>(Z,V);
				Vector T=CGAL::cross_product<K>(V,W);
				W=W/std::sqrt(W.squared_length());
				T=T/std::sqrt(T.squared_length());
				Point P(P1.x()*s+P2.x()*(1.-s),P1.y()*s+P2.y()*(1.-s),P1.z()*s+P2.z()*(1.-s));
				for(double u=-1;u<=1;u+=0.2)
				{
					for(double v=-1;v<=1;v+=0.2)
					{
						if(u*u+v*v<=1)
						{
							interior.push_back(P+Radius*W*u+Radius*T*v);
						}
					}
				}
			}
			itrC++;
			itrD++;
		}
		itrA++;
		itrB++;
	}
	
	// building their Delaunay tria*/ngulation.
	std::cout<<"start triangulation"<<endl;
	
	Delaunay T;
	//C2t3 ct(T);
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
	std::cout<<"T.number_of_cells:"<<T.number_of_cells()<<endl;
	std::cout<<"T.number_of_facets:"<<T.number_of_facets()<<endl;

	int max=*(std::max_element(cells.begin(),cells.end()));
	for(int i=0;i<=max;i++)
	{
		int n=std::count(cells.begin(),cells.end(),i);
		if(n!=0)
		std::cout<<i<<":"<<n<<endl;
	}
	/*Ž©‘OŽÀ‘•‚É‚·‚é*/
	NN=T.number_of_facets()/20;
	N=0;
	std::vector<Delaunay::Facet> facet_list;
	/*for(Delaunay::Finite_facets_iterator itr=T.finite_facets_begin();itr!=T.finite_facets_end();itr++,N++)
	{
		std::map<Delaunay::Cell_iterator,int>::iterator it_c1=index.find(itr->first);
		Delaunay::Cell_handle _cell=itr->first->neighbor(itr->second);
		std::map<Delaunay::Cell_iterator,int>::iterator it_c2=index.find(_cell);
		if(it_c2==index.end())
		{
			facet_list.push_back(*itr);
			std::cout<<"boundary"<<endl;
		}else
		{
			int d1=it_c1->second;
			int d2=it_c2->second;
			//if(cells[d1]==0&&cells[d2]<1)
			//{
				facet_list.push_back(*itr);
			//}
		}
		if(((int)N/NN)*NN==N)std::cout<<"*";
	}*/
	for(Delaunay::Cell_iterator itr=T.cells_begin();itr!=T.cells_end();itr++,N++)
	{
		std::map<Delaunay::Cell_iterator,int>::iterator it_c=index.find(itr);
		Delaunay::Cell_handle _cell0=itr->neighbor(0);
		Delaunay::Cell_handle _cell1=itr->neighbor(1);
		Delaunay::Cell_handle _cell2=itr->neighbor(2);
		Delaunay::Cell_handle _cell3=itr->neighbor(3);
		std::map<Delaunay::Cell_iterator,int>::iterator it_c0=index.find(_cell0);
		std::map<Delaunay::Cell_iterator,int>::iterator it_c1=index.find(_cell1);
		std::map<Delaunay::Cell_iterator,int>::iterator it_c2=index.find(_cell2);
		std::map<Delaunay::Cell_iterator,int>::iterator it_c3=index.find(_cell3);
		int D=6;
		if(cells[it_c->second]>=D/*&&cells[it_c0->second]<D*/)
		facet_list.push_back(Delaunay::Facet(itr,0));
		if(cells[it_c->second]>=D/*&&cells[it_c1->second]<D*/)
		facet_list.push_back(Delaunay::Facet(itr,1));
		if(cells[it_c->second]>=D/*&&cells[it_c2->second]<D*/)
		facet_list.push_back(Delaunay::Facet(itr,2));
		if(cells[it_c->second]>=D/*&&cells[it_c3->second]<D*/)
		facet_list.push_back(Delaunay::Facet(itr,3));
		if(((int)N/NN)*NN==N)std::cout<<"*";
	}
	std::cout<<endl;
	std::cout<<"complex created"<<endl;
	std::cout<<"T.number_of_facets:"<<T.number_of_facets()<<endl;


	//File write

	
	std::map<Delaunay::Vertex_handle,int> vIndex;
	filename="c:/out/"+NAME+".out";
	string filename2="c:/out/"+NAME+"_S.out";
	string filename3="c:/out/"+NAME+"_F.out";
	std::cout<<"start writing file"<<"["<<filename<<"]"<<endl;
	ofstream ofs(filename);
	ofstream ofs2(filename2);
	ofstream ofs3(filename3);
	for(vector<Point>::iterator itr=interior.begin();itr!=interior.end();itr++)
	{
		ofs3<<(*itr).x()<<" , "<<(*itr).y()<<" , "<<(*itr).z()<<endl;
	}
	N=0;
	for(vector<Delaunay::Facet>::iterator itr=facet_list.begin();itr!=facet_list.end();itr++)
	{
		for(int i=0;i<4;i++)
		{
			if(i!=itr->second)
			{
				Delaunay::Vertex_handle handle=itr->first->vertex(i);				
				std::map<Delaunay::Vertex_handle,int>::iterator pair=vIndex.find(handle);
				if(pair==vIndex.end())
				{
					Delaunay::Point center=itr->first->circumcenter();
					Delaunay::Point P=handle->point();
					Vector V=(P-center)*0.9;
					Delaunay::Point Q(center.x()+V.x(),center.y()+V.y(),center.z()+V.z());
					ofs<<Q.x()<<" , "<<Q.y()<<" , "<<Q.z()<<endl;
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
	ofs.close();
	ofs2.close();
	ofs3.close();
	

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
	for(vector<planes*>::iterator itr=planeTree.begin();itr!=planeTree.end();itr++)
	{
		delete(*itr);
	}
	//delete(cells);



  return 0;
}

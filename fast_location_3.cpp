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
	string NAME="donut";
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
	std::list<std::pair<Point,double>> interior;
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

    int NN=(data.size()/20);
	if(NN<1)NN=1;
	N=0;
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
		if(((int)N/NN)*NN==N)
		{
			std::cout<<"*";
		}
		N++;
	}
	std::cout<<endl;
	itrA=data.begin();
	itrB=planeTree.begin();
	std::cout<<"construct interior"<<endl;
	N=0;
    NN=(data.size()/20);
	if(NN<1)NN=1;

	while(itrA!=data.end())
	{
		double Radius;
		br* _branch;
		planes* _planes=*itrB;
		boost::tie(Radius,_branch)=*itrA;
		double dx=Radius/10.;
		br::iterator itrC=_branch->begin();
		planes::iterator itrD=_planes->begin();
		int N2=0;
		while(itrC!=_branch->end()-1)
		{
			//std::cout<<"N="<<N<<","<<"N2="<<N2<<endl;
			Vector V1=*itrD;
			Vector V2=*(itrD+1);
			Point P1=*itrC;
			Point P2=*(itrC+1);
			//•ªŠ„”Œˆ‚ß‚é
			double L=std::sqrt((P1-P2).squared_length());
			double DX=0.2;
			if(L<dx*5)
			{
				DX=0.2;
			}else
			{
				DX=dx/L;
			}
			for(double s=DX/2.;s<1.0;s+=DX)
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
				for(double u=-1;u<=1;u+=0.1)
				{
					for(double v=-1;v<=1;v+=0.1)
					{
						if(u*u+v*v<=1)
						{
							Point PP=P+Radius*W*u+Radius*T*v;
							interior.push_back(std::make_pair(PP,DX*L*0.1*Radius*0.1*Radius));
						}
					}
				}
			}
			itrC++;
			itrD++;
			N2++;
		}
		itrA++;
		itrB++;
		if(((int)N/NN)*NN==N)
		{
			std::cout<<"*";
		}
		N++;
	}
	std::cout<<endl;
	// building their Delaunay triangulation.
	std::cout<<"start triangulation"<<endl;
	
	Delaunay T;
	std::cout<<"exterior.size()"<<exterior.size()<<endl;
	std::cout<<"interior.size()"<<interior.size()<<endl;
	N=0;
	NN=(exterior.size()/20);
	if(NN<1)NN=1;
	for(int i=0;i<20;i++)
	{
		std::vector<Point>::iterator be=exterior.begin()+exterior.size()*i/20;
		std::vector<Point>::iterator en=exterior.begin()+exterior.size()*(i+1)/20;
		T.insert(be,en);
		std::cout<<"*";
	}
	std::cout<<endl;
	std::cout<<"end triangulation"<<endl;
	std::cout<<"T.number_of_vertices:"<<T.number_of_vertices()<<endl;
	std::cout<<"number_of_finite_cells:"<<T.number_of_finite_cells()<<endl;
	N=0;
	std::map<Delaunay::Cell_handle,int> index;
	//std::map<Delaunay::Cell_handle,double> volume;

	double* cells=new double[T.number_of_finite_cells()];
	//double* volume=new double[T.number_of_finite_cells()];
	double* density=new double[T.number_of_finite_cells()];
	//std::map<Delaunay::Cell_handle,double> density;
	std::cout<<"start cell search"<<endl;
	std::cout<<"create index"<<endl;

	for(Delaunay::Finite_cells_iterator itr=T.finite_cells_begin();itr!=T.finite_cells_end();itr++,N++)
	{
		index.insert(pair<const Delaunay::Finite_cells_iterator,int>(itr,N));
	}

	
	Cell_handle before;
	Cell_handle c;
	std::cout<<"initialize cells"<<endl;
	int S=T.number_of_finite_cells();
	double* ptr1=cells;
	double* ptr2=density;
	for(int i=0;i<S;i++)
	{
		*ptr1=0;
		*ptr2=0;
		ptr1++;
		ptr2++;
	}
	Locate_type lt;
	int li, lj;
	N=0;
	std::cout<<"cell size"<<endl;
	NN=(interior.size()/20);
	if(NN<1)NN=1;
	std::cout<<"start cell locate"<<endl;
	int S1=0,S2=0;
	for(list<std::pair<Point,double>>::iterator itr=interior.begin();itr!=interior.end();itr++,N++)
	{
		Point P=itr->first;
		if(N==0)
			c = T.locate(P, lt,li,lj);
		else
			c=T.locate(P,lt,li,lj,before);
		if(lt==Locate_type::CELL)
		{
				std::map<Delaunay::Cell_handle,int>::iterator it_c=index.find(c);
				int d=it_c->second;
				cells[d]+=itr->second;
				before=c;
				S1++;
		}
		if(((int)N/NN)*NN==N)
		{
			std::cout<<"*";
		}

	}
	std::cout<<endl;
	std::cout<<"normal="<<S1<<","<<"error="<<S2<<endl;
	std::cout<<"end cell search"<<endl;
	N=0;
	for(Delaunay::Finite_cells_iterator itr=T.finite_cells_begin();itr!=T.finite_cells_end();itr++,N++)
	{
		Point PA=itr->vertex(0)->point();
		Point PB=itr->vertex(1)->point();
		Point PC=itr->vertex(2)->point();
		Point PD=itr->vertex(3)->point();
		Vector G=PD-PA;
		Vector H=PD-PB;
		Vector I=PD-PC;
		double V=std::fabs(CGAL::cross_product(G,H)*I)/6.;
		density[N]=cells[N]/V;
	}

	std::cout<<"T.number_of_finite_cells:"<<T.number_of_finite_cells()<<endl;
	std::cout<<"T.number_of_finite_facets:"<<T.number_of_finite_facets()<<endl;

	NN=T.number_of_finite_facets()/20;
	if(NN<1)NN=1;
	N=0;
	std::vector<Delaunay::Facet> facet_list;
	for(Delaunay::Finite_cells_iterator itr=T.finite_cells_begin();itr!=T.finite_cells_end();itr++,N++)
	{
		double V=density[N];
		double D=0.2;
		for(int i=0;i<4;i++)
		{
			Delaunay::Cell_handle _neighbor=itr->neighbor(i);
			int c=index.find(_neighbor)->second;
			double Vi=density[c];
			std::map<Delaunay::Cell_handle,int>::iterator it_N=index.find(_neighbor);
			/*if(V>D)
			{
				std::cout<<V<<endl;
				facet_list.push_back(Delaunay::Facet(itr,i));
			}*/
			if(it_N==index.end())
			{
				if(V>D)
					facet_list.push_back(Delaunay::Facet(itr,i));
			}else
			{
				if(V>D&&Vi<D)
					facet_list.push_back(Delaunay::Facet(itr,i));
			
			}
		}
		if(((int)N/NN)*NN==N)std::cout<<"*";
	}
	std::cout<<endl;
	std::cout<<"complex created"<<endl;
	std::cout<<"T.number_of_finite_facets:"<<T.number_of_finite_facets()<<endl;

	
	//File write

	
	std::map<Delaunay::Vertex_handle,int> vIndex;
	filename="c:/out/"+NAME+".out";
	string filename2="c:/out/"+NAME+"_S.out";
	string filename3="c:/out/"+NAME+"_F.out";
	//string filename4="c:/out/"+NAME+"_D.out";
	std::cout<<"start writing file"<<"["<<filename<<"]"<<endl;
	ofstream ofs(filename);
	ofstream ofs2(filename2);
	ofstream ofs3(filename3);
	/*ofstream ofs4(filename4);
	for(list<std::pair<Point,double>>::iterator itr=interior.begin();itr!=interior.end();itr++)
	{
		ofs4<<(itr->first).x()<<" , "<<(itr->first).y()<<" , "<<(itr->first).z()<<endl;
	}*/
	N=0;
	for(Delaunay::Finite_cells_iterator itr=T.finite_cells_begin();itr!=T.finite_cells_end();itr++,N++)
	{
		Point PA=itr->vertex(0)->point();
		Point PB=itr->vertex(1)->point();
		Point PC=itr->vertex(2)->point();
		Point PD=itr->vertex(3)->point();
		double x=(PA.x()+PB.x()+PC.x()+PD.x())/4.;
		double y=(PA.y()+PB.y()+PC.y()+PD.y())/4.;
		double z=(PA.z()+PB.z()+PC.z()+PD.z())/4.;
		Point center(x,y,z);
		ofs3<<center.x()<<" , "<<center.y()<<" , "<<center.z()<<" , "<<density[N]<<endl;
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
					Delaunay::Point P=handle->point();
					ofs<<P.x()<<" , "<<P.y()<<" , "<<P.z()<<endl;
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
	//ofs4.close();
	
	
	std::cout << "Press Return To Exit...";
	std::cin.get();

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
	delete(cells);
	//delete(volume);
	delete(density);



  return 0;
}

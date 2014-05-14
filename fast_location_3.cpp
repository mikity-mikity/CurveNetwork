#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
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
//typedef CGAL::Triangulation_3<K>      Delaunay;
//typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned,K>    Vb;
//typedef CGAL::Triangulation_data_structure_3<Vb>                    Tds;
typedef CGAL::Delaunay_triangulation_3<K,CGAL::Fast_location> Delaunay;
typedef Delaunay::Point Point;
typedef CGAL::Vector_3<K> Vector;
typedef std::vector<Point> br;
typedef struct
{
	Vector majorRadius;
	Vector minorRadius;
	Vector Normal;
}eclipse;

typedef std::vector<eclipse> eclipses;
typedef boost::tuple<double,br*> R_br;
typedef Delaunay::Cell_handle Cell_handle;
typedef Delaunay::Vertex_handle Vertex_handle;
typedef Delaunay::Locate_type    Locate_type;
using namespace std;
int main()
{
	double pi=boost::math::constants::pi<double>();
	string NAME="test2";
	string filename="c:/out/"+NAME+".txt";
	std::cout<<"start reading file"<<"["<<filename<<"]"<<endl;

	ifstream ifs(filename);


	if(ifs.fail()){
		std::cout<<"File does not exist."<<endl;
		exit(0);
	}
	std::vector<R_br> data;
	std::vector<eclipses*> eclipseTree;
	br* branch;
	string line;
	unsigned N;
	double R;

	// read file
	double minR=10000;
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
			if(R<minR)minR=R;
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
	//basic resolution
	double baseRes=minR*2*pi/12.;

	std::vector<Point> exterior;
	std::list<std::pair<Point,double>> interior;
	int C=0;
	int num=0;
	
	//Construct eclipses at each node	
	for(vector<R_br>::iterator itr=data.begin();itr!=data.end();itr++,C++)
	{
		double Radius;
		br* _branch;
		boost::tie(Radius,_branch)=*itr;
		//planes *_planes=new planes();
		eclipses *_eclipses=new eclipses();
		eclipseTree.push_back(_eclipses);
		for(br::iterator anitr=_branch->begin();anitr!=_branch->end();anitr++)
		{
			eclipse newEclipse;
			if(anitr==_branch->begin())
			{
				Point P=*(anitr);
				Point Q=*(anitr+1);
				Vector V=(Q-P);
				V=V/std::sqrt(V.squared_length());
				Vector Z(0,0,1);
				double t=V*Z;
				if(std::fabs(t)>0.9)
				{
					Z=Vector(0,1,0);
				}
				Vector W=CGAL::cross_product(V,Z);
				Vector T=CGAL::cross_product(W,V);
				W=W/std::sqrt(W.squared_length());
				T=T/std::sqrt(T.squared_length());
				newEclipse.majorRadius=W;
				newEclipse.minorRadius=T;
				newEclipse.Normal=CGAL::cross_product(W,T);
			}else if(anitr==_branch->end()-1)
			{
				Point Q=*(anitr-1);
				Point R=*(anitr);
				Vector V=(R-Q);
				V=V/std::sqrt(V.squared_length());
				Vector Z(0,0,1);
				double t=V*Z;
				if(std::fabs(t)>0.9)
				{
					Z=Vector(0,1,0);
				}
				Vector W=CGAL::cross_product(V,Z);
				Vector T=CGAL::cross_product(W,V);
				W=W/std::sqrt(W.squared_length());
				T=T/std::sqrt(T.squared_length());
				newEclipse.majorRadius=W;
				newEclipse.minorRadius=T;
				newEclipse.Normal=CGAL::cross_product(W,T);
			}else
			{
				Point P=*(anitr-1);
				Point Q=*(anitr);
				Point R=*(anitr+1);
				Vector V1=(Q-P);
				Vector V2=(R-Q);
				V1=V1/std::sqrt(V1.squared_length());
				V2=V2/std::sqrt(V2.squared_length());
				if(V1*V2==1.0)//if parallel
				{
					Vector Z(0,0,1);
					double t=V1*Z;
					if(std::fabs(t)>0.9)
					{
						Z=Vector(0,1,0);
					}
					Vector W=CGAL::cross_product(V1,Z);
					Vector T=CGAL::cross_product(W,V1);
					W=W/std::sqrt(W.squared_length());
					T=T/std::sqrt(T.squared_length());
					newEclipse.majorRadius=W;
					newEclipse.minorRadius=T;
					newEclipse.Normal=CGAL::cross_product(W,T);
				}else
				{
					//angle between two lines
					double inner=V1*V2;
					double alpha=std::acos(inner);
					double theta=alpha/2.;

					Vector N=CGAL::cross_product(V1,V2);
					Vector E=CGAL::cross_product(V2,N);
					E=E/std::sqrt(E.squared_length());
					Vector R=(std::cos(theta)*E-std::sin(theta)*V2);
					newEclipse.majorRadius=R;
					newEclipse.minorRadius=N;
					newEclipse.Normal=CGAL::cross_product(R,N);
				}

			}
			_eclipses->push_back(newEclipse);
		}
	}
	vector<R_br>::iterator itrA=data.begin();
	vector<eclipses*>::iterator itrB=eclipseTree.begin();

	std::cout<<"construct exterior"<<endl;

    int NN=(data.size()/20);
	if(NN<1)NN=1;
	N=0;
	while(itrA!=data.end())
	{
		double Radius;
		br* _branch;
		eclipses* _eclipses=*itrB;
		boost::tie(Radius,_branch)=*itrA;
		br::iterator itrC=_branch->begin();
		eclipses::iterator itrD=_eclipses->begin();
		Vector beforeX(0,0,0);
		Vector beforeY(0,0,0);
		//Division number along diameter
		int RDIV=12;
		if(Radius*2.*pi/12.<baseRes)
		{
			RDIV=12;
		}else{
			RDIV=(int)(Radius*2.*pi/baseRes);
		}
		while(itrC!=_branch->end()-1)
		{
			//exterior.push_back(*itrC);
			//Division number along line
			Point P=*itrC;
			Point Q=*(itrC+1);
			Vector V=Q-P;
			double Length=std::sqrt(V.squared_length());
			V=V/Length;
			int DIV=1;//default division number
			if(Length<baseRes)
			{
				DIV=1;
			}else{
				DIV=(int)(Length/baseRes)+1;
			}
			//if in the first run, before is igonored
			//from the second run, before is used

			if(itrC==_branch->begin())
			{
				Vector Z(0,0,1);
				double t=V*Z;
				if(std::fabs(t)>0.9)
				{
					Z=Vector(0,1,0);
				}
				beforeX=CGAL::cross_product(V,Z);
				beforeY=CGAL::cross_product(V,beforeX);

			}else
			{				
				//project beforeX and beforeY to the plane at the node
				double cx=-(beforeX*itrD->Normal)/(V*itrD->Normal);
				double cy=-(beforeY*itrD->Normal)/(V*itrD->Normal);
				Vector tmpX=beforeX+cx*V;
				Vector tmpY=beforeY+cy*V;

				//Compute plane	
				Vector Z(0,0,1);
				double t=V*Z;
				if(std::fabs(t)>0.9)
				{
					Z=Vector(0,1,0);
				}
				Vector X=CGAL::cross_product(V,Z);
				Vector Y=CGAL::cross_product(V,beforeX);
				Vector N=CGAL::cross_product(X,Y);
				//project tmpX and tmpY to the plane conformed of X and Y
				beforeX=tmpX-(tmpX*N)*N;
				beforeY=CGAL::cross_product(V,beforeX);

			}
			for(int ss=0;ss<=DIV;ss++)
			{
				double s=((double)ss)/((double)DIV);
				exterior.push_back(Point((*itrC).x()*(1-s)+(*(itrC+1)).x()*s,(*itrC).y()*(1-s)+(*(itrC+1)).y()*s,(*itrC).z()*(1-s)+(*(itrC+1)).z()*s));
				for(int i=0;i<RDIV;i++)
				{
					double theta=(double)i/RDIV*2.*pi;

					Vector BI=0.8*Radius*(beforeX*std::cos(theta)+beforeY*std::sin(theta));
					Vector BE=1.05*Radius*(beforeX*std::cos(theta)+beforeY*std::sin(theta));

					//Project N
					double cI1=-(BI*itrD->Normal)/(V*itrD->Normal);
					double cI2=-(BI*(itrD+1)->Normal)/(V*(itrD+1)->Normal);
					double cE1=-(BE*itrD->Normal)/(V*itrD->Normal);
					double cE2=-(BE*(itrD+1)->Normal)/(V*(itrD+1)->Normal);
					Vector tmpI1=BI+cI1*V;
					Vector tmpI2=BI+cI2*V;
					Vector tmpE1=BE+cE1*V;
					Vector tmpE2=BE+cE2*V;
					Point DI1=(*itrC)+tmpI1;
					Point DI2=(*(itrC+1))+tmpI2;
					Point DE1=(*itrC)+tmpE1;
					Point DE2=(*(itrC+1))+tmpE2;
					if(itrC==_branch->end()-2)
					{
						Point DI(DI2.x()*s+DI1.x()*(1-s),DI2.y()*s+DI1.y()*(1-s),DI2.z()*s+DI1.z()*(1-s));
						Point DE(DE2.x()*s+DE1.x()*(1-s),DE2.y()*s+DE1.y()*(1-s),DE2.z()*s+DE1.z()*(1-s));
						exterior.push_back(DI);
						exterior.push_back(DE);
					}else
					{
						Point DI(DI2.x()*s+DI1.x()*(1-s),DI2.y()*s+DI1.y()*(1-s),DI2.z()*s+DI1.z()*(1-s));
						Point DE(DE2.x()*s+DE1.x()*(1-s),DE2.y()*s+DE1.y()*(1-s),DE2.z()*s+DE1.z()*(1-s));
						exterior.push_back(DI);
						exterior.push_back(DE);
					}
				}
			}
			itrC++;
			itrD++;
		}
		//exterior.push_back(*(_branch->end()-1));

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
	std::cout<<"isValid"<<T.is_valid()<<endl;
	std::cout<<"T.number_of_vertices:"<<T.number_of_vertices()<<endl;
	std::cout<<"number_of_finite_cells:"<<T.number_of_finite_cells()<<endl;
	N=0;
	std::map<Delaunay::Cell_handle,int> index;
	//std::map<Delaunay::Cell_handle,double> volume;

	int* cells=new int[T.number_of_finite_cells()];
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
	int* ptr1=cells;
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
	
	std::cout<<"start cell locate"<<endl;
	

	N=0;
	itrA=data.begin();
	itrB=eclipseTree.begin();
	while(itrA!=data.end())
	{
		int YDIV=8;
		double Radius;
		br* _branch;
		eclipses* _eclipses=*itrB;
		boost::tie(Radius,_branch)=*itrA;
		br::iterator itrC=_branch->begin();
		eclipses::iterator itrD=_eclipses->begin();
		Vector beforeX(0,0,0);
		Vector beforeY(0,0,0);
		//Division number along diameter
		int RDIV=12;
		if(Radius*2.*pi/12.<baseRes)
		{
			RDIV=12;
		}else{
			RDIV=(int)(Radius*2.*pi/baseRes);
		}

		//Generate base particles
		vector<Vector> vectors;
		vector<Point> particles;
		double alpha=pi/((double)RDIV);
		double R2=Radius*std::cos(alpha);
		for (int i=0;i<RDIV;i++)
		{
			double theta=2.*pi*((double)i)/((double)RDIV);
			Vector vector(R2*std::cos(theta),R2*std::sin(theta),0);
			vectors.push_back(vector);
		}
		double edge=2*std::sin(alpha)*Radius;
		double dx=edge/YDIV;
		for(double i=-Radius;i<Radius;i+=dx)
		{
			for(double j=-Radius;j<Radius;j+=dx)
			{
				Vector V(i,j,0);
				bool flag=true;
				for(auto T=vectors.begin();T!=vectors.end();T++)
				{
					if ((*T)*V<R2*R2)
					{
						continue;
					}else
					{
						flag=false;
						break;
					}
            
				}
				if(flag)
					particles.push_back(Point(i,j,0));
			}
		}

		while(itrC!=_branch->end()-1)
		{
			//Division number along line
			Point P=*itrC;
			Point Q=*(itrC+1);
			Vector V=Q-P;
			double Length=std::sqrt(V.squared_length());
			V=V/Length;
			int DIV=1;//default division number
			if(Length<baseRes)
			{
				DIV=1;
			}else{
				DIV=(int)(Length/baseRes)+1;
			}
			DIV*=5;

			for(double ss=0.5;ss<DIV;ss++)
			{					
				for(auto particle=particles.begin();particle!=particles.end();particle++)
				{
						N++;
				}
			}
			itrC++;
			itrD++;
		}
		itrA++;
		itrB++;
	}
	std::cout<<"interior.size():"<<N<<endl;
	int size=N;
	int next=1000000;
	N=0;
	itrA=data.begin();
	itrB=eclipseTree.begin();
	while(itrA!=data.end())
	{
		int YDIV=8;
		double Radius;
		br* _branch;
		eclipses* _eclipses=*itrB;
		boost::tie(Radius,_branch)=*itrA;
		br::iterator itrC=_branch->begin();
		eclipses::iterator itrD=_eclipses->begin();
		Vector beforeX(0,0,0);
		Vector beforeY(0,0,0);
		//Division number along diameter
		int RDIV=12;
		if(Radius*2.*pi/12.<baseRes)
		{
			RDIV=12;
		}else{
			RDIV=(int)(Radius*2.*pi/baseRes);
		}

		//Generate base particles
		vector<Vector> vectors;
		vector<Point> particles;
		double alpha=pi/((double)RDIV);
		double R2=Radius*std::cos(alpha);
		for (int i=0;i<RDIV;i++)
		{
			double theta=2.*pi*((double)i)/((double)RDIV);
			Vector vector(R2*std::cos(theta),R2*std::sin(theta),0);
			vectors.push_back(vector);
		}
		double edge=2*std::sin(alpha)*Radius;
		double dx=edge/YDIV;
		for(double i=-Radius;i<Radius;i+=dx)
		{
			for(double j=-Radius;j<Radius;j+=dx)
			{
				Vector V(i,j,0);
				bool flag=true;
				for(auto T=vectors.begin();T!=vectors.end();T++)
				{
					if ((*T)*V<R2*R2)
					{
						continue;
					}else
					{
						flag=false;
						break;
					}
            
				}
				if(flag)
					particles.push_back(Point(i,j,0));
			}
		}

		while(itrC!=_branch->end()-1)
		{
			exterior.push_back(*itrC);
			//Division number along line
			Point P=*itrC;
			Point Q=*(itrC+1);
			Vector V=Q-P;
			double Length=std::sqrt(V.squared_length());
			V=V/Length;
			int DIV=1;//default division number
			if(Length<baseRes)
			{
				DIV=1;
			}else{
				DIV=(int)(Length/baseRes)+1;
			}
			DIV*=5;
			//if in the first run, before is igonored
			//from the second run, before is used

			if(itrC==_branch->begin())
			{
				Vector Z(0,0,1);
				double t=V*Z;
				if(std::fabs(t)>0.9)
				{
					Z=Vector(0,1,0);
				}
				beforeX=CGAL::cross_product(V,Z);
				beforeY=CGAL::cross_product(V,beforeX);

			}else
			{				
				//project beforeX and beforeY to the plane at the node
				double cx=-(beforeX*itrD->Normal)/(V*itrD->Normal);
				double cy=-(beforeY*itrD->Normal)/(V*itrD->Normal);
				Vector tmpX=beforeX+cx*V;
				Vector tmpY=beforeY+cy*V;

				//Compute plane	
				Vector Z(0,0,1);
				double t=V*Z;
				if(std::fabs(t)>0.9)
				{
					Z=Vector(0,1,0);
				}
				Vector X=CGAL::cross_product(V,Z);
				Vector Y=CGAL::cross_product(V,beforeX);
				Vector N=CGAL::cross_product(X,Y);
				//project tmpX and tmpY to the plane conformed of X and Y
				beforeX=tmpX-(tmpX*N)*N;
				beforeY=CGAL::cross_product(V,beforeX);

			}
			for(double ss=0.5;ss<DIV;ss++)
			{
				double s=((double)ss)/((double)DIV);
				for(auto particle=particles.begin();particle!=particles.end();particle++)
				{

					Vector B=beforeX*(*particle).x()+beforeY*(*particle).y();

					//Project N
					double c1=-(B*itrD->Normal)/(V*itrD->Normal);
					double c2=-(B*(itrD+1)->Normal)/(V*(itrD+1)->Normal);
					Vector tmp1=B+c1*V;
					Vector tmp2=B+c2*V;
					Point D1=(*itrC)+tmp1;
					Point D2=(*(itrC+1))+tmp2;
					Point D(D2.x()*s+D1.x()*(1-s),D2.y()*s+D1.y()*(1-s),D2.z()*s+D1.z()*(1-s));
					if(N==0)
						c = T.locate(D, lt,li,lj);
					else
						c=T.locate(D,lt,li,lj,before);
						
					if(lt==Locate_type::CELL)
					{
							std::map<Delaunay::Cell_handle,int>::iterator it_c=index.find(c);
							int d=it_c->second;
							cells[d]++;
							before=c;
					}
					if(N>=next)
					{
						std::cout<<N<<"/"<<size<<endl;
						next+=1000000;
					}
					N++;
				}
			}

			itrC++;
			itrD++;
		}
		itrA++;
		itrB++;
	}
	std::cout<<endl;
	
	
	
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
		density[N]=V;//cells[N]/V;
	}
	
	std::cout<<"T.number_of_finite_cells:"<<T.number_of_finite_cells()<<endl;
	std::cout<<"T.number_of_finite_facets:"<<T.number_of_finite_facets()<<endl;
	bool* bool_list=new bool[T.number_of_finite_cells()];
	bool* bool_list2=new bool[T.number_of_finite_cells()];
	N=0;
	double D=0.4;
	int count=0;
	for(Delaunay::Finite_cells_iterator itr=T.finite_cells_begin();itr!=T.finite_cells_end();itr++,N++)
	{
		bool_list[N]=true;
		if(cells[N]<1)bool_list[N]=false;
		bool_list2[N]=true;
	}
	for(int i=0;i<200;i++)
	{
		N=0;
		for(Delaunay::Finite_cells_iterator itr=T.finite_cells_begin();itr!=T.finite_cells_end();itr++,N++)
		{
			if(bool_list2[N]==true)
			{
				bool flag=false;
				for(int i=0;i<4;i++)
				{
					Delaunay::Cell_handle _neighbor=itr->neighbor(i);
					std::map<Delaunay::Cell_handle,int>::iterator it_N=index.find(_neighbor);
					if(it_N==index.end())
					{
						flag=true;break;
					}else
					{
						if(bool_list2[N]&&(!bool_list2[it_N->second]))
						{
							flag=true;break;
						}
					}
				}
				if(flag)
				{
					if(bool_list[N]==false)bool_list2[N]=false;
				}
			}
		}
		std::cout<<"*";
	}
	std::cout<<endl;
	N=0;


	NN=T.number_of_finite_facets()/20;
	if(NN<1)NN=1;
	N=0;
	std::vector<Delaunay::Facet> facet_list;
	for(Delaunay::Finite_cells_iterator itr=T.finite_cells_begin();itr!=T.finite_cells_end();itr++,N++)
	{
		for(int i=0;i<4;i++)
		{
			Delaunay::Cell_handle _neighbor=itr->neighbor(i);
			std::map<Delaunay::Cell_handle,int>::iterator it_N=index.find(_neighbor);
			if(it_N==index.end())
			{
				if(bool_list2[N])
					facet_list.push_back(Delaunay::Facet(itr,i));
			}else
			{
				if(bool_list2[N]&&(!bool_list2[it_N->second]))
				{
					//if(density[it_N->second]>0.00001)
					facet_list.push_back(Delaunay::Facet(itr,i));
				}
			}
		}
		if(((int)N/NN)*NN==N)std::cout<<"*";
	}
	std::cout<<endl;
	std::cout<<"complex created"<<endl;
	std::cout<<"T.number_of_finite_facets:"<<T.number_of_finite_facets()<<endl;
	delete(bool_list);
	delete(bool_list2);
	
	//File write

	
	std::map<Delaunay::Vertex_handle,int> vIndex;
	filename="c:/out/"+NAME+".out";
	string filename2="c:/out/"+NAME+"_S.out";
	//string filename3="c:/out/"+NAME+"_F.out";
	//string filename4="c:/out/"+NAME+"_D.out";
	std::cout<<"start writing file"<<"["<<filename<<"]"<<endl;
	ofstream ofs(filename);
	ofstream ofs2(filename2);
	//ofstream ofs3(filename3);
	/*ofstream ofs4(filename4);
	for(list<std::pair<Point,double>>::iterator itr=interior.begin();itr!=interior.end();itr++)
	{
		ofs4<<(itr->first).x()<<" , "<<(itr->first).y()<<" , "<<(itr->first).z()<<endl;
	}*/
	N=0;
	/*for(Delaunay::Finite_cells_iterator itr=T.finite_cells_begin();itr!=T.finite_cells_end();itr++,N++)
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
	}*/
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
	//ofs3.close();
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
	for(vector<eclipses*>::iterator itr=eclipseTree.begin();itr!=eclipseTree.end();itr++)
	{
		delete(*itr);
	}
	delete(cells);
	//delete(volume);
	delete(density);



  return 0;
}

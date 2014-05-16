#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/IO/Color.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_vertex_base_with_info_3<CGAL::Color, K> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb>                    Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>                      Delaunay;

typedef Delaunay::Point   Point;
using namespace std;
int main()
{
	Delaunay T;

	T.insert(Point(0,0,0));
	T.insert(Point(1,0,0));
	T.insert(Point(0,1,0));
	T.insert(Point(0,0,1));
	T.insert(Point(1./4.,1./4.,1./4.));
	std::map<Delaunay::Cell_iterator,unsigned> index;
	int N=0;
	for(auto itr=T.cells_begin();itr!=T.cells_end();itr++,N++)
	{
		index.insert(std::pair<const Delaunay::Cell_iterator,unsigned>(itr,N));
	}
	for(auto itr=T.cells_begin();itr!=T.cells_end();itr++,N++)
	{
		std::cout<<"N:"<<N<<endl;
		std::cout<<itr->vertex(0)->point().x()<<","<<itr->vertex(0)->point().y()<<","<<itr->vertex(0)->point().z()<<endl;
		std::cout<<itr->vertex(1)->point().x()<<","<<itr->vertex(1)->point().y()<<","<<itr->vertex(1)->point().z()<<endl;
		std::cout<<itr->vertex(2)->point().x()<<","<<itr->vertex(2)->point().y()<<","<<itr->vertex(2)->point().z()<<endl;
		std::cout<<itr->vertex(3)->point().x()<<","<<itr->vertex(3)->point().y()<<","<<itr->vertex(3)->point().z()<<endl;
		Delaunay::Cell_handle cell0=itr->neighbor(0);
		Delaunay::Cell_handle cell1=itr->neighbor(1);
		Delaunay::Cell_handle cell2=itr->neighbor(2);
		Delaunay::Cell_handle cell3=itr->neighbor(3);
		int i1=0,i2=0,i3=0,i4=0;
		if(index.find(cell0)!=index.end())
		{
			i1=index.find(cell0)->second;
		}else i1=1000;
		if(index.find(cell1)!=index.end())
		{
			i2=index.find(cell1)->second;
		}else i2=1000;
		if(index.find(cell2)!=index.end())
		{
			i3=index.find(cell2)->second;
		}else i3=1000;
		if(index.find(cell3)!=index.end())
		{
			i4=index.find(cell3)->second;
		}else i4=1000;
		std::cout<<i1<<","<<i2<<","<<i3<<","<<i4<<endl;
	}
	// Set the color of finite vertices of degree 6 to red.
	std::cout<<"T.number_of_vertices"<<T.number_of_vertices()<<endl;
	std::cout<<"T.number_of_cells"<<T.number_of_cells()<<endl;
	std::cout<<"T.number_of_finite_cells"<<T.number_of_finite_cells()<<endl;
	std::cout<<"T.number_of_facets"<<T.number_of_facets()<<endl;
	std::cout<<"T.number_of_finite_facets"<<T.number_of_finite_facets()<<endl;
	std::cout << "Press Return To Exit...";
	std::cin.get();

  return 0;
}

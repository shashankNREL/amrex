/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */


#ifndef _STLEXPLORER_H_
#define _STLEXPLORER_H_

#include "AMReX_RealVect.H"
#include "AMReX_IntVect.H"
#include "AMReX_Box.H"
#include "AMReX_STLMesh.H"
#include "AMReX_STLBox.H"
#include "AMReX_CellEdge.H"
#include "AMReX_KDTree.H"

using std::shared_ptr;
namespace amrex
{

/*
 * Class to explore an STL mesh
 * has member functions that
 * 1) build a K-D tree of the vertices
 * 2) search for nearest triangles given
 *  a) a point (return 1 triangle)
 *  b) a cell/volume (return all triangles with some part inside the volume)
 * 3) find intersection of a triangle with a single edge
 * 4) find intersection(s) of the whole mesh and an edge
 */

  class STLExplorer
  {
  public:

    /// Constructor - just need a mesh
    STLExplorer(shared_ptr<STLMesh> a_stlmesh);

    /*
   /// Copy constructor
   STLExplorer(const STLExplorer& a_inputExplorer);
    */

    /// Destructor
    ~STLExplorer();

    /// builds cellToTriangles - connectivity between box and stlmesh
    void Explore(shared_ptr<STLBox> a_sb);

    void Explore(const Box&           a_region,
                 const Box          & a_domain,
                 const RealVect&      a_origin,
                 const RealVect&      a_dx);


    /// return the point of intersection between a cellEdge and the mesh
    /// and whether the two nodes are inside or outside
    /// a_intersectPt is bogus if both nodes are inside or both outside
    void GetCellEdgeIntersection(const CellEdge& a_celledge,
                                 RealVect&       a_intersectPt,
                                 bool&           a_isNode0Inside,
                                 bool&           a_isNode1Inside);

    /// return whether or not a point is inside or outside the domain
    void GetPointInOut(const IntVect&  a_point,
                       bool&           a_inout);

    // get methods for temporary variables
    void GetVertMap(Vector<IntVect>** a_vertmap);
    void GetTriMap(Vector<Vector<IntVect> >** a_trimap);
    void GetSTLBox(shared_ptr<STLBox>& a_sb);

  protected:

    // actual data (shared)
    shared_ptr<STLMesh> m_msh;   // pointer to mesh
    shared_ptr<STLBox>  m_sb;    // pointer to mesh<->box data
    KDTree*                m_ptree; // pointer to the KDTree of nodes

    bool m_freestlbox; // true if STLExplorer should free the stlbox when destroyed
    bool m_printdebug; // if true, will print (lots of) debug info

    // (temporary data holder) vertex index -> cell (cell that a vertex is in)
    // only used in FindCellsOnEdges (to fill in cells between vertices of mesh)
    // NOTE: this structure may contain IntVects that are invalid (outside of m_box)
    Vector<IntVect> m_vertmap;

    // (temporary data holder) triangle index -> cell (Vector of cells in a given triangle)
    // only used in FindCellsInTriangles (to fill in intervening cells within triangles)
    // NOTE: this structure may contain IntVects that are invalid (outside of m_box)
    Vector<Vector<IntVect> > m_trimap;

    // overall function to explore a box
    void DoExplore();

    // sets map<IntVect,TriInCell> cellmap (adds on cells as we find them)
    // and  Vector<Intvect> vertmap
    void FindCellsOnVertices();

    // goes along triangle edges and finds cells (sets left & right triangles)
    // sets cellmap (inserts cells) and trimap
    void FindCellsOnEdges();

    // fills in any cells completely contained by the triangle
    // sets cellmap (inserts cells) and trimap
    void FindCellsInTriangles();

    // removes all cells that lie outside of the Box from consideration
    // modifies cellmap, trimap, vertmap (call before setting edgemap and nodemap)
    void RemoveCellsOutsideDomain();

    // stores the set of cell edges on the boundary
    // builds edgemap and nodemap
    void FindCellEdgesOnBondary();

    // builds a KDTree of the nodes that we know are inside/outside
    void BuildKDTree();

    // returns whether node0 or node1 of the cell edge is inside the domain
    bool WhichNodeIsInside(const CellEdge& celledge,
                           const int&      triangle);

    // returns whether node0 and node1 
    void FindEdgeInOut(const CellEdge& celledge,
                       bool&           isNode0Inside,
                       bool&           isNode1Inside);

    // same as above, but uses a KDTree instead of exhaustive search
    void FindEdgeInOutWithKDTree(const CellEdge& celledge,
                                 bool&           isNode0Inside,
                                 bool&           isNode1Inside);

    // given a triangle and an edge, find the intersection
    RealVect FindPlaneLineIntersection(const CellEdge& celledge,
                                       const int&      triangle);

    // return true if point lies within the triangle (up to the tolerance mesh.tol??)
    bool IsPointInTriangle(const RealVect& point,
                           const int&      triangle);

    // return true if the point lies on the edge (up to the tolerance mesh.tol)
    bool IsPointOnCellEdge(const RealVect& point,
                           const CellEdge& celledge);

    // note, i is an input and is incremented in FillInCellLine
    void FillInCellLine(vector<IntVect>& cells,
                        int&             i,
                        const int&       itri,
                        const int&       idir0,
                        const int&       idir1);

  protected:
    void RecursiveKDTreeInsert(vector<pair<RealVect,pair<IntVect,bool>*> > &allNodes,
                               const int                                   &nstart,
                               const int                                   &nend,
                               const int                                   &depth);

  private:
    STLExplorer()
      {
      }

    void operator=(const STLExplorer& a_inputReader)
      {
      }
  };
}

#endif


#include <vector>
class SparseMatrix{
  private:
    int rowSize, colSize;
    std::vector<int> ptr, idx; 
    std::vector<double> a;
    void setStorage(const std::vector<double>& s);
    void setPtr(const std::vector<int>& ptr);
    void setIdx(const std::vector<int>& idx);
  public:
    SparseMatrix();
    ~SparseMatrix();
    void setSize(int, int);
    void setCCS(const std::vector<std::vector<double> >& mat);
    void setParam(const std::vector<double>& stora, const std::vector<int>& ptr, const std::vector<int>& idx);
    std::vector<std::vector<double> > transSparseXmat(const std::vector<std::vector<double> >& );
    std::vector<std::vector<double> > matXsparse(const std::vector<std::vector<double> >& );
};

#include <vector>
class SparseMatrix{
  private:
    int rowSize, colSize;
    std::vector<int> ptr, idx; 
    std::vector<double> a;
    std::vector<double> at;
  public:
    SparseMatrix();
    ~SparseMatrix();
    void setSize(int, int);
    void setStorage(const std::vector<double>& s);
    void setParam(const std::vector<int>& ptr, const std::vector<int>& idx);
    std::vector<std::vector<double> > transSparseXmat(const std::vector<std::vector<double> >& );
    std::vector<std::vector<double> > matXsparse(const std::vector<std::vector<double> >& );
};

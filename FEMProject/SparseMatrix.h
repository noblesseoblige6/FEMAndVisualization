#include <vector>
class SparseMatrix{
  private:
    int rowSize, colSize;
    std::vector<int> rowPtr, rowIdx; 
    std::vector<int> colPtr, colIdx; 
    std::vector<double> a;
    std::vector<double> at;
  public:
    SparseMatrix();
    ~SparseMatrix();
    void setSize(int, int);
    void setStorage(const std::vector<double>& s);
    void setRowParam(const std::vector<int>& ptr, const std::vector<int>& idx);
    void setColParam(const std::vector<int>& ptr, const std::vector<int>& idx);
    void transposeParam();
    std::vector<std::vector<double> > sparseXmat(const std::vector<std::vector<double> >& );
    std::vector<std::vector<double> > matXsparse(const std::vector<std::vector<double> >& );
};

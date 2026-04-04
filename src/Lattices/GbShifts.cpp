//
// Created by Nikhil Chandra Admal on 2/4/24.
//
#include "../../include/Lattices/GbShifts.h"
#include "../../include/Utilities/randomInteger.h"

namespace oILAB {
template <int dim>
GbShifts<dim>::GbShifts(const Gb<dim> &gb,
                        const ReciprocalLatticeVector<dim>& axis,
                        const std::vector<LatticeVector<dim>>& gbCslVectors,
                        const double& tMax,
                        const double& sPerpMax) try :
/* init */ gb(gb)
/* init */,axis(axis)
/* init */,gbCslVectors(gbCslVectors)
/* init */,tShiftPairs(getbShiftPairs(gb, gbCslVectors, tMax, sPerpMax))
{
    std::cout << "--------------------GBShifts class construction "
                 "---------------------------"
              << std::endl;
    std::cout << "GB CSL vectors = " << std::endl;
    for (const auto& elem : gbCslVectors)
        std::cout << elem.cartesian().transpose() << std::endl;
    std::cout << std::endl;

    double latticeConstant= gb.bc.A.latticeBasis.col(0).norm();
    std::cout << "Maximum translation = "
              << tMax * latticeConstant << std::endl;
    std::cout << std::endl;
    std::cout << "Maximum shift offset = "
              << sPerpMax * latticeConstant << std::endl;
    std::cout << std::endl;

    VectorDimD normal;
    if (dim == 3)
        normal = gbCslVectors[0].cross(gbCslVectors[1]).cartesian().normalized();
    else
        normal = gbCslVectors[0].cross().cartesian().normalized();

    Eigen::IOFormat fmt(12, 0, " ", " ", "", "", "", "");
    std::cout << std::fixed << std::setprecision(6);


    std::cout << "Exploring the following translation-shift pairs:" << std::endl;
    for (const auto &[t, s] : tShiftPairs) {
        std::cout << "t = " << t.cartesian().transpose().format(fmt);
        std::cout << "; s = " << s.transpose().format(fmt) << std::endl;
        assert(t.cartesian().norm() <= tMax*latticeConstant && "Translation exceeds the max translation.\n");
        assert(s.dot(normal) <= sPerpMax*latticeConstant/2.0 && "Shift exceeds the max offset along the GB normal.\n");
    }
    std::cout << "----------------------------" << std::endl;
    std::cout << std::endl;
}
catch(std::runtime_error& e)
{
    std::cout << e.what() << std::endl;
    throw(std::runtime_error("GB construction failed. "));
}

template<int dim>
std::vector<std::pair<LatticeVector<dim>, typename GbShifts<dim>::VectorDimD>> GbShifts<dim>::getbShiftPairs(const Gb<dim>& gb,
                                                                                                             const std::vector<LatticeVector<dim>>& gbCslVectors,
                                                                                                             const double& tMax,
                                                                                                             const double& sPerpMax)
{
    std::vector<std::pair<LatticeVector<dim>, VectorDimD>> output;

    // ensure the input CSL vectors describe the GB.
    assert(gbCslVectors.size()==dim-1);

    // form the CSL cell for modulo operations to identify the shift vectors
    // - the cell is spanned by the input GB CSL vectors and an out-of-plane vector scaled
    //   by a factor determined by sPerpMax
    auto nC= gb.bc.getReciprocalLatticeDirectionInC(gb.nB.reciprocalLatticeVector());
    auto gbPlaneParallelCslBasis= gb.bc.csl.planeParallelLatticeBasis(nC,true);
    std::vector<LatticeVector<dim>> cslSubLatticeVectors;
    double latticeConstant= gb.bc.A.latticeBasis.col(0).norm();
    int factor= floor(sPerpMax*latticeConstant/nC.planeSpacing() +FLT_EPSILON);
    factor= (factor>0 ? factor : 1);
    cslSubLatticeVectors.push_back(factor*gbPlaneParallelCslBasis[0].latticeVector());
    cslSubLatticeVectors.push_back(gbCslVectors[0]);
    cslSubLatticeVectors.push_back(gbCslVectors[1]);
    auto cslPoints= gb.bc.csl.box(cslSubLatticeVectors,"cslSubLattice.txt");

    // Collect all DSCL translations within a parallelopiped of size tMax
    //  - orient this parallelopiped along the tilt axis, period vector directions,
    //    and the third direction is close to being parallel to the GB normal
    std::vector<LatticeVector<dim>> orthogonalCslLatticeVectors;
    for (int i=0; i<dim; ++i)
        orthogonalCslLatticeVectors.push_back(cslSubLatticeVectors[i]);
    gb.bc.updateBoxVectors(orthogonalCslLatticeVectors,0.8);
    std::vector<LatticeVector<dim>> latticeVectorsDscl;
    for (int i=0; i<dim; ++i)
        latticeVectorsDscl.push_back(gb.bc.getLatticeDirectionInD(orthogonalCslLatticeVectors[i]).latticeVector());
    for(int i=0; i<dim; ++i)
    {
        int factor= floor(tMax*latticeConstant/latticeVectorsDscl[i].cartesian().norm()+FLT_EPSILON);
        factor= (factor>0 ? factor : 1);
        latticeVectorsDscl[i]= factor*latticeVectorsDscl[i];
    }
    auto allTranslations= gb.bc.dscl.box(latticeVectorsDscl,"translations.txt");

    VectorDimD shiftT, shiftC;
    shiftT << -0.5, -0.5, -0.5;
    shiftC << -0.5, -FLT_EPSILON, -FLT_EPSILON;
    for(auto& translation: allTranslations) {
        LatticeVector<dim>::modulo(translation, latticeVectorsDscl, shiftT);
        if (translation.cartesian().norm() > tMax*gb.bc.A.latticeBasis.col(0).norm())
            continue;
        auto cslShift = LatticeVector<dim>((gb.bc.LambdaA * translation).eval(), gb.bc.dscl);
        for(const auto& cslPoint : cslPoints) {
            VectorDimD cslShiftCentered = cslShift.cartesian() + cslPoint.cartesian() - translation.cartesian() / 2;
            LatticeVector<dim>::modulo(cslShiftCentered, cslSubLatticeVectors, shiftC);
            if (abs(cslShiftCentered.dot(nC.cartesian().normalized())) <= sPerpMax*latticeConstant/2.0)
                output.push_back(std::make_pair(translation, cslShiftCentered));
        }
    }
    return output;
}
/*
template<int dim>
std::vector<LatticeVector<dim>> GbShifts<dim>::getGbCslVectors(const Gb<dim>& gb, const ReciprocalLatticeVector<dim>& axis)
{
    std::vector<LatticeVector<dim>> output;
    LatticeVector<dim> axisA(gb.bc.A.latticeDirection(axis.cartesian()).latticeVector());
    LatticeVector<dim> axisC(gb.bc.getLatticeDirectionInC(axisA).latticeVector());
    for(int i=0; i<dim-1; ++i) {
        if (i==0) output.push_back(gb.getPeriodVector(axis));
        if (i==1) output.push_back(axisC);
    }
    return output;
}
 */

//template class GbShifts<2>;
template class GbShifts<3>;

} // namespace oILAB

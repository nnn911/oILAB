/* This file is part of gbLAB.
 *
 * gbLAB is distributed without any warranty under the MIT License.
 */


#ifndef gbLAB_BiCrystal_h_
#define gbLAB_BiCrystal_h_

#include "../Math/RationalMatrix.h"
#include "../Math/SmithDecomposition.h"
#include "LatticeModule.h"
#include <array>
#include <map>

namespace oILAB {

    /*! \brief Character of a grain boundary generated relative to a misorientation axis:
     * Tilt (axis lies in the GB plane), Twist (axis is normal to the GB plane), or
     * Mixed (axis is at a general angle to the GB plane).
     */
    enum class GBCharacter { Tilt, Twist, Mixed };

    /*!Class template that computes the coincident-site-lattice (CSL) of two
     * parent lattices using the Smith Normal Form [1].
     *
     * [1] Coincidence Lattices and Associated Shear Transformations
     */
    template <int dim>
    class BiCrystal : public RationalMatrix<dim>
    /*             */,public SmithDecomposition<dim>
    {
        typedef typename LatticeCore<dim>::IntScalarType IntScalarType;
        typedef typename LatticeCore<dim>::VectorDimD VectorDimD;
        typedef typename LatticeCore<dim>::MatrixDimD MatrixDimD;
        typedef typename LatticeCore<dim>::VectorDimI VectorDimI;
        typedef typename LatticeCore<dim>::MatrixDimI MatrixDimI;

        
        static MatrixDimI getM(const RationalMatrix<dim>& rm, const SmithDecomposition<dim>& sd);
        static MatrixDimI getN(const RationalMatrix<dim>& rm, const SmithDecomposition<dim>& sd);
        static MatrixDimI getLambdaA(const MatrixDimI& M, const MatrixDimI& N);
        static MatrixDimI getLambdaB(const MatrixDimI& M, const MatrixDimI& N);
        static MatrixDimD getCSLBasis(const Lattice<dim>& A,
                                      const Lattice<dim>& B,
                                      const SmithDecomposition<dim>& sd,
                                      const MatrixDimI& M,
                                      const MatrixDimI& N);
        static MatrixDimD getDSCLBasis(const Lattice<dim>& A,
                                       const Lattice<dim>& B,
                                       const SmithDecomposition<dim>& sd,
                                       const MatrixDimI& M,
                                       const MatrixDimI& N);

    public:
        template<IntScalarType keyScale>
        class GBKey
        {
        public:
            GBKey(const Gb<dim>& gb) : _referenceVector(gb.nA.cartesian().normalized()) {}

            [[nodiscard]] IntScalarType operator()(const Gb<dim>& gb) const noexcept
            {
                const double cosAngle = gb.nA.cartesian().normalized().dot(_referenceVector);
                return (std::clamp(cosAngle, -1.0, 1.0) + 1.0) * keyScale;
            }

            [[nodiscard]] constexpr static IntScalarType minKey() noexcept { return 0; };
            [[nodiscard]] constexpr static IntScalarType maxKey() noexcept { return 2 * keyScale; };
            [[nodiscard]] constexpr static IntScalarType numKeys() noexcept { return maxKey() + 1; };

        private:
            ReciprocalLatticeVector<dim>::VectorDimD _referenceVector;
        };

    public:
        const Lattice<dim>& A;
        const Lattice<dim>& B;

        /*! \brief Integer matrix that connects the bases \f$ \textbf A^{\|} \f$ and \f$ \textbf C^{\|} \f$ of
         * lattices \f$\mathcal A\f$ and the CSL \f$\mathcal C\f$, respectively: \f$ \textbf C^{\|} = \textbf A^{\|} \textbf M\f$.
         */
        const MatrixDimI M;

        /*! \brief Integer matrix that connects the bases \f$ \textbf B^{\|} \f$ and \f$ \textbf C^{\|} \f$ of
         * lattices \f$\mathcal B\f$ and the CSL \f$\mathcal C\f$, respectively: \f$ \textbf C^{\|} = \textbf B^{\|} \textbf N\f$.
         */
        const MatrixDimI N;

        /*! \brief Signed ratio of the unit cell volume of \f$\mathcal C\f$ to that of \f$\mathcal A\f$.
         *  \f$ \Sigma_{\mathcal A} = \det(\textbf M)\f$.
         */
        const IntScalarType sigmaA;

        /*! \brief Signed ratio of the unit cell volume of \f$\mathcal C\f$ to that of \f$\mathcal B\f$:
         *  \f$ \Sigma_{\mathcal B} = \det(\textbf N)\f$.
         */
        const IntScalarType sigmaB;

        /*! \brief \f$ \Sigma = |\Sigma_{\mathcal A}|\f$ if \f$ \Sigma_{\mathcal A} = \Sigma_{\mathcal B}\f$, else
         *  \f$ \Sigma = 0 \f$
         */
        const int sigma;

    private:
        /*! \brief Lattice \f$\mathcal C\f$ with basis \f$\textbf C^\|\f$
         */
        const Lattice<dim> cslp;

        /*! \brief Lattice \f$\mathcal D\f$ with basis \f$\textbf D^\|\f$
         */
        const Lattice<dim> dsclp;

        /*! \brief Maps integer coordinates w.r.t \f$\textbf D\f$ basis to
         * integer coordinates w.r.t parallel basis \f$\textbf D^\|\f$
         */
        const MatrixDimI dscl2dsclp;

        /*! \brief Maps integer coordinates w.r.t \f$\textbf C\f$ basis to
         * integer coordinates w.r.t parallel basis \f$\textbf C^\|\f$
         */
        const MatrixDimI csl2cslp;

    public:
        /*! \brief CSL lattice \f$\mathcal C\f$.
         */
        const Lattice<dim> csl;

        /*! \brief DCSL lattice \f$\mathcal D\f$.
         */
        const Lattice<dim> dscl;

        /*! \brief Lattice \f$\mathcal A\f$ with basis \f$\textbf A^\|\f$
         */
        const Lattice<dim> Ap;

        /*! \brief Lattice \f$\mathcal B\f$ with basis \f$\textbf B^\|\f$
         */
        const Lattice<dim> Bp;

        /*! \brief Shift tensor \f$\Lambda_{\mathcal A}:\mathbb Z_{\mathcal D} \to \mathbb Z_{\mathcal D}\f$ describes the
         * shift in the CSL when lattice \f$\mathcal A\f$ is shifted by a DSCL vector.
         */
        const MatrixDimI LambdaA;

        /*! \brief Shift tensor \f$\Lambda_{\mathcal B}:\mathbb Z_{\mathcal D} \to \mathbb Z_{\mathcal D}\f$ describes the
         * shift in the CSL when lattice \f$\mathcal B\f$ is shifted by a DSCL vector.
         */
        const MatrixDimI LambdaB;


        LatticeVector<dim> shiftTensorA(const LatticeVector<dim>& d) const;
        LatticeVector<dim> shiftTensorB(const LatticeVector<dim>& d) const;

        /**********************************************************************/
        /*! \brief Constructs a bicrystal from two lattices \f$\mathcal A \f$ and \f$\mathcal B \f$ by computing the
         *  parallel bases \f$\textbf A^\|, \textbf B^\|, \textbf C^\|, \textbf D^\|\f$ for lattices
         *  \f$\mathcal A, \mathcal B, \mathcal C, \mathcal D\f$, respectively.
         *  If the flag useRLLL is .true., then the bases of CSL and DSCL are reduced using the LLL algorithm.
         *
         * \param[in] Lattices  \f$\mathcal A \f$ and \f$\mathcal B \f$, and useRLLL flag
         * \returns   A bicrystal object
         * */
        BiCrystal(const Lattice<dim>& A,
                  const Lattice<dim>& B,
                  const bool& useRLLL=false);

        /*!
         * Outputs lattice vector in lattice \f$\mathcal A\f$ that is equal to the inputted vector \f$\textbf v\f$
         * that belongs to \f$\mathcal A\f$ or \f$\mathcal C\f$
         * @param v - lattice vector
         * @return LatticeVector in \f$\mathcal A\f$
         */
        LatticeVector<dim> getLatticeVectorInA(const LatticeVector<dim>& v) const;
        /*!
         * Outputs lattice vector in lattice \f$\mathcal B\f$ that is equal to the inputted vector \f$\textbf v\f$
         * that belongs to \f$\mathcal B\f$ or \f$\mathcal C\f$
         * @param v - lattice vector
         * @return LatticeVector in \f$\mathcal B\f$
         */
        LatticeVector<dim> getLatticeVectorInB(const LatticeVector<dim>& v) const;
        /*!
         * Outputs lattice vector in lattice \f$\mathcal D\f$ that is equal to the inputted vector \f$\textbf v\f$
         * that belongs to \f$\mathcal A\f$, \f$\mathcal B\f$, \f$\mathcal C\f$, or \f$\mathcal D\f$
         * @param v - lattice vector
         * @return LatticeVector in \f$\mathcal D\f$
         */
        LatticeVector<dim> getLatticeVectorInD(const LatticeVector<dim>& v) const;
        /*!
         * Outputs lattice direction in the CSL \f$\mathcal C\f$ that is parallel to the inputted vector \f$\textbf v\f$
         * that belongs to one of the four lattices, \f$\mathcal A\f$, \f$\mathcal B\f$, \f$\mathcal C\f$, or \f$\mathcal D\f$,
         * @param v - lattice vector
         * @return LatticeDirection in \f$\mathcal C\f$
         */
        LatticeDirection<dim> getLatticeDirectionInC(const LatticeVector<dim>& v) const;
        /*!
         * Outputs lattice direction in the DSCL \f$\mathcal D\f$ that is parallel to the inputted vector \f$\textbf v\f$
         * that belongs to one of the four lattices, \f$\mathcal A\f$, \f$\mathcal B\f$, \f$\mathcal C\f$, or \f$\mathcal D\f$,
         * @param v - lattice vector
         * @return LatticeDirection in \f$\mathcal D\f$
         */
        LatticeDirection<dim> getLatticeDirectionInD(const LatticeVector<dim>& v) const;

        /*!
         * Outputs reciprocal lattice direction in the dual lattice \f$\mathcal A^*\f$ that is parallel to the inputted
         * reciprocal vector \f$\textbf v\f$ that belongs to one of the four dual lattices, \f$\mathcal A^*\f$, \f$\mathcal B^*\f$,
         * \f$\mathcal C^*\f$, or \f$\mathcal D^*\f$,
         * @param v - lattice vector
         * @return ReciprocalLatticeDirection in \f$\mathcal A^*\f$
         */
        ReciprocalLatticeDirection<dim> getReciprocalLatticeDirectionInA(const ReciprocalLatticeVector<dim>& v) const;
        /*!
         * Outputs reciprocal lattice direction in the dual lattice \f$\mathcal B^*\f$ that is parallel to the inputted
         * reciprocal vector \f$\textbf v\f$ that belongs to one of the four dual lattices, \f$\mathcal A^*\f$, \f$\mathcal B^*\f$,
         * \f$\mathcal C^*\f$, or \f$\mathcal D^*\f$,
         * @param v - lattice vector
         * @return ReciprocalLatticeDirection in \f$\mathcal B^*\f$
         */
        ReciprocalLatticeDirection<dim> getReciprocalLatticeDirectionInB(const ReciprocalLatticeVector<dim>& v) const;
        /*!
         * Outputs reciprocal lattice direction in the dual lattice \f$\mathcal C^*\f$ that is parallel to the inputted
         * reciprocal vector \f$\textbf v\f$ that belongs to one of the four dual lattices, \f$\mathcal A^*\f$, \f$\mathcal B^*\f$,
         * \f$\mathcal C^*\f$, or \f$\mathcal D^*\f$,
         * @param v - lattice vector
         * @return ReciprocalLatticeDirection in \f$\mathcal C^*\f$
         */
        ReciprocalLatticeDirection<dim> getReciprocalLatticeDirectionInC(const ReciprocalLatticeVector<dim>& v) const;
        /*!
         * Outputs reciprocal lattice direction in the dual lattice \f$\mathcal D^*\f$ that is parallel to the inputted
         * reciprocal vector \f$\textbf v\f$ that belongs to one of the four dual lattices, \f$\mathcal A^*\f$, \f$\mathcal B^*\f$,
         * \f$\mathcal C^*\f$, or \f$\mathcal D^*\f$,
         * @param v - lattice vector
         * @return ReciprocalLatticeDirection in \f$\mathcal D^*\f$
         */
        ReciprocalLatticeDirection<dim> getReciprocalLatticeDirectionInD(const ReciprocalLatticeVector<dim>& v) const;

        /*!
         * \brief Given an axis \f$\textbf d\f$, that belongs to lattices \f$\mathcal A\f$ or \f$\mathcal B\f$, this
         * function generates a set of tilt, twist, or mixed GBs, depending on \p character. CURRENTLY ONLY WORDS FOR DIMENSION 3
         * @param d - LatticeDirection that describes the misorientation axis
         * @param div - parameter to span the GBs (ignored when character == GBCharacter::Twist)
         * @param character - GBCharacter::Tilt (default), GBCharacter::Twist, or GBCharacter::Mixed
         * @return A data structure that stores GBs sorted (primarily) in increasing order of their inclination angle.
         */
        std::map<IntScalarType,Gb<dim>>
        generateGrainBoundaries(const LatticeDirection<dim>& d, int div=30, GBCharacter character = GBCharacter::Tilt) const requires (dim==2 || dim==3);

        /*!
         * \brief Given an axis \f$\textbf d\f$, that belongs to lattices \f$\mathcal A\f$ or \f$\mathcal B\f$, this
         * function generates a set of tilt, twist, or mixed GBs, depending on \p character. CURRENTLY ONLY WORDS FOR DIMENSION 3
         * @param d - LatticeDirection that describes the misorientation axis
         * @param div - parameter to span the GBs (ignored when character == GBCharacter::Twist)
         * @param character - GBCharacter::Tilt (default), GBCharacter::Twist, or GBCharacter::Mixed
         * @param callback - callback function that is called for each GB
         */
        template<typename Callback>
        void generateGrainBoundaries(const LatticeDirection<dim>& d, int div = 30, GBCharacter character = GBCharacter::Tilt, Callback&& callback = {}) const
            requires(dim == 2 || dim == 3);

        /*!
         * \brief Refines \p boxVectors in place so that the box they define is as orthogonal as
         * possible, trading off against \p orthogonality. Every vector remains a CSL lattice
         * vector; vectors may grow (become a larger, equivalent linear combination) but the
         * plane/line each one spans together with the others is unchanged, so \p boxVectors
         * continues to describe the same GB.
         *
         * For \p dim==3, \p boxVectors[1] and \p boxVectors[2] are assumed to already span the
         * GB plane. \p boxVectors[2] is left untouched; \p boxVectors[1] is first replaced by
         * the in-plane CSL vector that is as orthogonal as possible to \p boxVectors[2] (found
         * by recursing into this same search one dimension down, on the trivial \f$\Sigma=1\f$
         * bicrystal of the plane's own 2D CSL).
         *
         * In both \p dim==2 and \p dim==3, \p boxVectors[0] is then replaced by the CSL vector
         * that is as orthogonal as possible to the plane/line spanned by the remaining box
         * vectors, by searching increasingly larger candidates along that plane's/line's own
         * reciprocal direction.
         *
         * \p orthogonality controls how long each of the two searches above keeps looking for a
         * more orthogonal candidate before settling: 0 accepts the first (smallest) candidate,
         * 1 searches the full range the plane's/line's periodicity admits. Larger values can
         * substantially increase the size of the resulting box.
         *
         * @param boxVectors - \p dim linearly independent CSL lattice vectors, modified in place.
         * @param orthogonality - a value in the interval \f$[0,1]\f$.
         */
        void
        updateBoxVectors(std::array<LatticeVector<dim>, dim>& boxVectors,
                              const double& orthogonality=0.0) const requires (dim==2 || dim==3);

        /*! This function outputs/prints a 2D bicrystal (two lattices that form the GB and
         * the CSL) bounded by a box defined using
         * two input box vectors. The box vectors have to be linearly independent lattice
         * vectors. The function optimizes boxVectors[0]
         * to make the box as orthogonal as possible depending on the \p orthogonality parameter.
         *
         *
         * @param boxVectors two linearly independent lattice vectors.
         * @param orthogonality (double) a value in the interval \f$[0,1]\f$.
         * @param filename (optional) name of the output file
         * @param orient (optional) While printing to a file, orient the system such that one of the box sides
         * is along the global x axis. This flag does not
         * influence the returning configuration, only the configuration printed to the file.
         * @return lattice points of the bicrystal (along with the CSL) bounded by the box (std::vector<LatticeVector<2>>).
         */
        std::vector<LatticeVector<dim>>
        box(const std::array<LatticeVector<dim>, dim>& boxVectors,
            const int& dsclFactor,
            std::string filename= "",
            bool orient=false) const requires (dim==2 || dim==3);
    };
    
    
} // end namespace

#include "BiCrystalImplementation.h"

#endif


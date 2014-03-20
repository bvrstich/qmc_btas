//nog enkele definities:
#include "SpinQuantum.h"
namespace btas { typedef SpinQuantum Quantum; };

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

#include <btas/common/blas_cxx_interface.h>
#include <btas/common/TVector.h>

#include <btas/DENSE/TArray.h>
#include <btas/DENSE/TConj.h>
#include <btas/DENSE/TCONTRACT.h>

#include <btas/DENSE/DArray.h>
#include <btas/DENSE/ZArray.h>

#include <btas/SPARSE/STConj.h>

#include <btas/QSPARSE/QSTArray.h>
#include <btas/QSPARSE/QSTLAPACK.h>

#include <btas/QSPARSE/QSZArray.h>

#include <MPSblas.h>

#include "Random.h"
#include "SpinHamiltonian.h"
#include "Coupling.h"
#include "Trotter.h"
#include "Walker.h"
#include "Tools.h"

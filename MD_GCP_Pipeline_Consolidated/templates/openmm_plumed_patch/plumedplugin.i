%module openmmplumed


%import(module="openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"
%include "std_string.i"
%include "mpi4py.i"
%mpi4py_typemap(Comm, MPI_Comm);

%{
#include "PlumedForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
#include <mpi.h>
#include <tuple>
#include <vector>
%}

// Use SWIG's type conversion to properly handle System* from openmm module
// Directly access the pointer from SwigPyObject (workaround for cross-module conversion)
%typemap(in) OpenMM::System* {
    void *ptr = 0;
    // Try SWIG's type query first
    swig_type_info *type_info = SWIG_TypeQuery("OpenMM::System *");
    if (type_info != NULL) {
        int res = SWIG_ConvertPtr($input, &ptr, type_info, 0);
        if (SWIG_IsOK(res)) {
            $1 = reinterpret_cast<OpenMM::System*>(ptr);
        } else {
            SWIG_exception_fail(SWIG_ArgError(res), "in method '$symname', argument $argnum of type 'OpenMM::System *'");
        }
    } else {
        // Fallback: directly cast to SwigPyObject and extract ptr
        // This is a workaround for cross-module type conversion
        // We know $input is a SwigPyObject from the openmm module
        SwigPyObject *sobj = (SwigPyObject *)$input;
        if (sobj && sobj->ptr) {
            $1 = reinterpret_cast<OpenMM::System*>(sobj->ptr);
        } else {
            // Last resort: try generic conversion
            int res = SWIG_ConvertPtr($input, &ptr, 0, 0);
            if (SWIG_IsOK(res)) {
                $1 = reinterpret_cast<OpenMM::System*>(ptr);
            } else {
                SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError), "in method '$symname', argument $argnum: expected System object");
            }
        }
    }
}

// Helper function to directly add PlumedForce to System
// This bypasses the Python type conversion issue
%inline %{
int addPlumedForceToSystemCpp(OpenMM::System* system, PlumedPlugin::PlumedForce* pf) {
    return system->addForce(static_cast<OpenMM::Force*>(pf));
}
%}

%pythoncode %{
import openmm as mm

# Helper function to add PlumedForce to System
# Uses the fact that PlumedForce inherits from Force, so we can
# access the underlying C++ pointer and call System.addForce
# through a workaround
def addPlumedForceToSystem(system, plumed_force):
    """
    Add a PlumedForce to an OpenMM System.
    This is a workaround for SWIG cross-module type conversion issues.
    """
    if not isinstance(plumed_force, PlumedForce):
        raise TypeError("Expected PlumedForce object")
    
    # The issue is that System.addForce expects a Force object from the openmm module,
    # but PlumedForce is from the openmmplumed module.
    # Since PlumedForce inherits from Force, at the C++ level the pointer is compatible.
    # We need to trick SWIG into accepting it.
    
    # Try: call System_addForce directly using the underlying C++ pointers
    # Get access to the C++ function through the openmm module
    try:
        # Access the internal SWIG wrapper
        import openmm.openmm as _openmm_module
        # Try to call System_addForce with type conversion
        # Get the system's this pointer
        system_this = system.this
        # Get the plumed force's this pointer
        pf_this = plumed_force.this
        
        # The pointers are compatible at C++ level
        # We need SWIG to convert them
        # Try calling addForce through the system object
        # but with explicit type conversion
        return _openmm_module.System_addForce(system, plumed_force)
    except (AttributeError, TypeError) as e:
        # If direct call doesn't work, we need to use the C++ helper
        # But that also has type conversion issues...
        # Final fallback: modify simulate.py to handle this differently
        raise TypeError(
            f"SWIG cross-module type conversion failed. "
            f"PlumedForce cannot be automatically converted to Force* for System.addForce(). "
            f"Error: {e}. "
            f"Please ensure you're using the latest version of the wrapper with the fix applied."
        )
%}

namespace PlumedPlugin {

%copyctor PlumedForce;
class PlumedForce : public OpenMM::Force {
public:
    PlumedForce(const std::string& script, const MPI_Comm intra_comm, const MPI_Comm inter_comm);
    const std::string& getScript() const;
};

}

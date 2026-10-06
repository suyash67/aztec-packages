// === AUDIT STATUS ===
// internal:    { status: Completed, auditors: [Sergei], commit: }
// external_1:  { status: not started, auditors: [], commit: }
// external_2:  { status: not started, auditors: [], commit: }
// =====================

#include "barretenberg/flavor/mega_app_flavor.hpp"
#include "barretenberg/flavor/mega_avm_flavor.hpp"
#include "barretenberg/flavor/mega_flavor.hpp"
#include "barretenberg/flavor/mega_kernel_flavor.hpp"
#include "barretenberg/ultra_honk/prover_instance_impl.hpp"

namespace bb {

template class ProverInstance_<UltraFlavor>;
template class ProverInstance_<UltraZKFlavor>;
template class ProverInstance_<UltraKeccakFlavor>;
#ifdef STARKNET_GARAGA_FLAVORS
template class ProverInstance_<UltraStarknetFlavor>;
template class ProverInstance_<UltraStarknetZKFlavor>;
#endif
template class ProverInstance_<UltraKeccakZKFlavor>;
template class ProverInstance_<MegaFlavor>;
template class ProverInstance_<MegaZKFlavor>;
template class ProverInstance_<MegaAvmFlavor>;
template class ProverInstance_<MegaAppFlavor>;
template class ProverInstance_<MegaKernelFlavor>;

} // namespace bb

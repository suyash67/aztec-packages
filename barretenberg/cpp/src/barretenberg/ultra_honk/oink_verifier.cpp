// === AUDIT STATUS ===
// internal:    { status: Completed, auditors: [Sergei], commit: }
// external_1:  { status: not started, auditors: [], commit: }
// external_2:  { status: not started, auditors: [], commit: }
// =====================

#include "barretenberg/ext/starknet/flavor/ultra_starknet_flavor.hpp"
#include "barretenberg/ext/starknet/flavor/ultra_starknet_zk_flavor.hpp"
#include "barretenberg/flavor/mega_app_flavor.hpp"
#include "barretenberg/flavor/mega_app_recursive_flavor.hpp"
#include "barretenberg/flavor/mega_avm_recursive_flavor.hpp"
#include "barretenberg/flavor/mega_flavor.hpp"
#include "barretenberg/flavor/mega_kernel_flavor.hpp"
#include "barretenberg/flavor/mega_kernel_recursive_flavor.hpp"
#include "barretenberg/flavor/mega_recursive_flavor.hpp"
#include "barretenberg/flavor/mega_zk_recursive_flavor.hpp"
#include "barretenberg/flavor/ultra_keccak_zk_flavor.hpp"
#include "barretenberg/flavor/ultra_zk_recursive_flavor.hpp"
#include "barretenberg/ultra_honk/oink_verifier_impl.hpp"

namespace bb {

template class OinkVerifier<UltraFlavor>;
template class OinkVerifier<UltraZKFlavor>;
template class OinkVerifier<UltraKeccakFlavor>;
#ifdef STARKNET_GARAGA_FLAVORS
template class OinkVerifier<UltraStarknetFlavor>;
template class OinkVerifier<UltraStarknetZKFlavor>;
#endif
template class OinkVerifier<UltraKeccakZKFlavor>;
template class OinkVerifier<MegaFlavor>;
template class OinkVerifier<MegaZKFlavor>;

// Recursive flavor instantiations
template class OinkVerifier<UltraRecursiveFlavor_<UltraCircuitBuilder>>;
template class OinkVerifier<UltraRecursiveFlavor_<MegaCircuitBuilder>>;
template class OinkVerifier<MegaRecursiveFlavor_<UltraCircuitBuilder>>;
template class OinkVerifier<MegaRecursiveFlavor_<MegaCircuitBuilder>>;
template class OinkVerifier<MegaZKRecursiveFlavor_<MegaCircuitBuilder>>;
template class OinkVerifier<MegaZKRecursiveFlavor_<UltraCircuitBuilder>>;
template class OinkVerifier<MegaAvmRecursiveFlavor_<UltraCircuitBuilder>>;
template class OinkVerifier<UltraZKRecursiveFlavor_<UltraCircuitBuilder>>;
template class OinkVerifier<UltraZKRecursiveFlavor_<MegaCircuitBuilder>>;
template class OinkVerifier<MegaAppFlavor>;
template class OinkVerifier<MegaKernelFlavor>;
template class OinkVerifier<MegaAppRecursiveFlavor>;
template class OinkVerifier<MegaKernelRecursiveFlavor>;

} // namespace bb

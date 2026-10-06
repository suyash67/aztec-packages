// === AUDIT STATUS ===
// internal:    { status: Planned, auditors: [Khashayar], commit: }
// external_1:  { status: not started, auditors: [], commit: }
// external_2:  { status: not started, auditors: [], commit: }
// =====================

#include "barretenberg/commitment_schemes/small_subgroup_ipa/small_subgroup_ipa_impl.hpp"
#include "barretenberg/commitment_schemes/utils/test_settings.hpp"
#include "barretenberg/constants.hpp"
#include "barretenberg/ecc/curves/bn254/bn254.hpp"
#include "barretenberg/eccvm/eccvm_flavor.hpp"
#include "barretenberg/eccvm/eccvm_short_monomial_flavor.hpp"
#include "barretenberg/eccvm/eccvm_translation_data.hpp"
#include "barretenberg/ext/starknet/flavor/ultra_starknet_zk_flavor.hpp"
#include "barretenberg/flavor/mega_zk_flavor.hpp"
#include "barretenberg/flavor/ultra_keccak_zk_flavor.hpp"
#include "barretenberg/flavor/ultra_zk_flavor.hpp"
#include "barretenberg/polynomials/polynomial.hpp"
#include "barretenberg/polynomials/univariate.hpp"
#include "barretenberg/stdlib/primitives/curves/grumpkin.hpp"
#include "barretenberg/sumcheck/zk_sumcheck_data.hpp"
#include "barretenberg/translator_vm/translator_flavor.hpp"
#include <array>
#include <type_traits>
#include <vector>

namespace bb {

// Instantiate with ZK Flavors
template class SmallSubgroupIPAProver<ECCVMFlavor>;
template class SmallSubgroupIPAProver<ECCVMShortMonomialFlavor>;
template class SmallSubgroupIPAProver<TranslatorFlavor>;
template class SmallSubgroupIPAProver<TranslatorShortMonomialFlavor>;
template class SmallSubgroupIPAProver<MegaZKFlavor>;
template class SmallSubgroupIPAProver<UltraZKFlavor>;
template class SmallSubgroupIPAProver<UltraKeccakZKFlavor>;
#ifdef STARKNET_GARAGA_FLAVORS
template class SmallSubgroupIPAProver<UltraStarknetZKFlavor>;
#endif

// Instantiations used in tests
template class SmallSubgroupIPAProver<BN254Settings>;
template class SmallSubgroupIPAProver<GrumpkinSettings>;

} // namespace bb

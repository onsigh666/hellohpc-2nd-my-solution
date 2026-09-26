#include "maimoe/kernel_api.hpp"

#include "MaiMoeEngine.hpp"

#include <cstddef>
#include <cstdint>
#include <span>
#include <string_view>

namespace maimoe::kernel {

// The engine is deterministic and maimoe_check is the authority on correctness,
// so this wrapper does no re-scanning or re-hashing of its own. Writing through
// process_chart_into also avoids the per-chart heap allocations and the extra
// copy that engine::process_chart + flatten would cost.
void process_chart(std::uint64_t chart_id, std::string_view chart_text,
                   std::span<std::uint8_t, kEncodedChartBytes> output) {
    engine::EncodedChartLayout layout;
    engine::process_chart_into(chart_id, chart_text, output, layout);
}

}

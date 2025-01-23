// problems/ProblemRegistry.cpp
#include "ProblemRegistry.hpp"
#include "Problem.hpp" // Include definitions of Problem and its derivatives
#include "SodShockTube.hpp"
#include "GravityCollapse.hpp"

namespace agoge
{
    namespace problems
    {

        std::unique_ptr<Problem> createProblem(const std::string &name)
        {
            if (name == "sod")
            {
                return std::make_unique<SodShockTube>();
            }
            else if (name == "collapse")
            {
                return std::make_unique<GravityCollapse>();
            }
            // Add more problem types as needed

            // If unknown:
            return nullptr;
        }

    } // namespace problems
} // namespace agoge

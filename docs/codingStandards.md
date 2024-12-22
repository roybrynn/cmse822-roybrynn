---

layout: default

title: Coding Standards

---

# Coding Standards

This page contains the coding standards for this course. The purpose of these standards is to create code and documentation that is well-structured, easy to understand, and maintainable. These standards are based on those for CMSE 201 but are updated to reflect the needs of this class and modern C++ practices.

--------------------------------------------------------------------------------

## Overriding principle for this coding standard

There is a saying popular among professional software developers, which we hope you will take to heart in this course and throughout your career:

"Always code as if the person who ends up maintaining your code will be an axe-wielding maniac who knows where you live."

In other words, write and document your code so that it is easy to read and understand. This is critical in the era of open-source scientific software—it helps others use (and reuse) your code and helps you do the same when revisiting your code months later.

Consistency and clarity are the most important goals of these standards. These guidelines are not meant to be rigid rules but rather a framework for writing clear, maintainable code.

--------------------------------------------------------------------------------

## Source code format

- **Line Length:** Source code lines should be no more than 79 characters long. Break long lines into shorter lines for readability.

- **Single Statement per Line:** Only one statement should be included per line of source code—do not put multiple statements in a single line.

- **Whitespace:** Use blank spaces and blank lines to enhance readability. Avoid excessive whitespace, which can decrease readability. For indentation, use **four spaces** (do NOT use tabs). Configure your editor (e.g., Visual Studio Code, CLion, emacs, vi, or nano) to use spaces instead of tabs.

- **Control Structures:** Avoid deep nesting of control structures. If clarity is impaired by too many levels of nesting, consider refactoring the deeply nested statements into separate functions.

- **Function Length:** Keep functions concise and focused on a single task. If a function becomes excessively long (e.g., hundreds of lines), break it into smaller functions to improve clarity and maintainability.

--------------------------------------------------------------------------------

## Naming conventions

- **General Principle:** Names should clearly describe their purpose and not be reused for multiple purposes. 

- **Variables:** Use "lower_snake_case" style for variable names. For example:

  ```cpp
  int car_count = 0; // Keeps track of cars entering the intersection
  ```

- **Constants:** Use symbolic constants rather than embedding arbitrary values directly in your code. Use all caps for constants:

  ```cpp
  const double BOLTZMANN_CONSTANT = 1.38e-23; // Boltzmann's constant in Joules/Kelvin
  double energy = BOLTZMANN_CONSTANT * temperature;
  ```

- **Functions:** Use the same "lower_snake_case" style for function names:

  ```cpp
  double calculate_energy(double temperature);
  ```

- **Classes:** Use "PascalCase" for class names:

  ```cpp
  class RocketShipEngine;
  ```

- **Enumerations:** Use "PascalCase" for enumeration names, and "ALL_CAPS" for enumeration values:

  ```cpp
  enum class EngineState { IDLE, RUNNING, STOPPED };
  ```

--------------------------------------------------------------------------------

## Descriptive comments

- **File Header Comments:** Every source code file must include a comment block at the top describing the file's purpose and contents:

  ```cpp
  /*
  File: conversion_utils.cpp
  Purpose: Contains utility functions for unit conversion.

  This file includes the following routines:
  - cm_to_parsecs: Converts centimeters to parsecs
  - oz_to_solar_masses: Converts ounces to solar masses
  - sec_to_millennia: Converts seconds to millennia
  */
  ```

- **Function Documentation:** Each function must include a block comment (or Doxygen-style comment) describing its purpose, parameters, and return values. For example:

  ```cpp
  /**
   * @brief Counts the number of cars passing through an intersection.
   *
   * @param intersection_number Integer describing the intersection number.
   * @param start_time Floating-point time to start counting.
   * @param stop_time Floating-point time to stop counting.
   * @return Integer representing the number of cars passing through the intersection.
   */
  int count_cars(int intersection_number, double start_time, double stop_time);
  ```

- **Inline Comments:** Use inline comments sparingly to explain complex or non-obvious code. Place comments above the line of code or block they describe, rather than at the end of the line.

--------------------------------------------------------------------------------

## Modern C++ Practices

- **Smart Pointers:** Use `std::unique_ptr` or `std::shared_ptr` instead of raw pointers for managing dynamically allocated memory.

  ```cpp
  std::unique_ptr<RocketShipEngine> engine = std::make_unique<RocketShipEngine>();
  ```

- **Standard Library Containers:** Prefer `std::vector`, `std::array`, and other STL containers over raw arrays.

  ```cpp
  std::vector<int> car_counts;
  ```

- **Range-based Loops:** Use range-based `for` loops when iterating over containers for better readability:

  ```cpp
  for (const auto& car : car_counts) {
      std::cout << car << std::endl;
  }
  ```

- **Error Handling:** Use exceptions for error handling instead of error codes, and use `try`/`catch` blocks judiciously to handle expected errors.

  ```cpp
  try {
      engine.start();
  } catch (const EngineException& e) {
      std::cerr << e.what() << std::endl;
  }
  ```

- **Default Initialization:** Initialize variables and objects immediately to avoid undefined behavior:

  ```cpp
  int car_count = 0;
  ```

- **Avoid Macros:** Use `constexpr` or `inline` functions instead of macros for constant values and inline logic.

  ```cpp
  constexpr double pi = 3.141592653589793;
  inline double area_of_circle(double radius) { return pi * radius * radius; }
  ```

--------------------------------------------------------------------------------

## References

[1] <http://www.cse.msu.edu/~cse231/General/coding.standard.html>

[2] <http://www.python.org/dev/peps/pep-0008/>

[3] "Introduction to Scientific Programming" by Victor Eijkhout


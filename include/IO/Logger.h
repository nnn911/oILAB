/*
 * This file is part of oILAB.
 * Copyright 2026 - 2026, Nikhil Chandra Admal and the oILAB contributors.
 * SPDX-License-Identifier: MIT
 */

#ifndef OILAB_LOGGER_H
#define OILAB_LOGGER_H

#include "TerminalColors.h"
#include <iostream>
#include <sstream>
#include <string>

namespace oILAB {

#ifdef OILAB_LOGGING_DISABLED

// Compile-time logging disabled: all Logger calls become no-ops with zero
// overhead.
struct NullStream {
  template <typename T> NullStream &operator<<(T &&) { return *this; }
};

enum class LogLevel { Debug = 0, Info = 1, Warn = 2, Error = 3 };

struct Logger {
  static LogLevel minLevel;
  static bool enabled;
  static NullStream debug() { return {}; }
  static NullStream info() { return {}; }
  static NullStream warn() { return {}; }
  static NullStream error() { return {}; }
};

inline LogLevel Logger::minLevel = LogLevel::Info;
inline bool Logger::enabled = true;

#else

enum class LogLevel { Debug = 0, Info = 1, Warn = 2, Error = 3 };

class LogStream {
  LogLevel level;
  bool active;
  std::ostringstream oss;

  static const char *prefix(LogLevel l) {
    switch (l) {
    case LogLevel::Debug:
      return "[debug] ";
    case LogLevel::Info:
      return "[info]  ";
    case LogLevel::Warn:
      return "[warn]  ";
    case LogLevel::Error:
      return "[error] ";
    }
    return "";
  }

  static const std::string &color(LogLevel l) {
    switch (l) {
    case LogLevel::Debug:
      return defaultColor;
    case LogLevel::Info:
      return cyanColor;
    case LogLevel::Warn:
      return yellowColor;
    case LogLevel::Error:
      return redBoldColor;
    }
    return defaultColor;
  }

public:
  LogStream(LogLevel l, bool a) : level(l), active(a) {}

  LogStream(LogStream &&o) noexcept : level(o.level), active(o.active) {
    oss << o.oss.str();
    o.active = false;
  }

  LogStream(const LogStream &) = delete;
  LogStream &operator=(const LogStream &) = delete;

  ~LogStream() {
    if (!active)
      return;
    std::string msg = oss.str();
    while (!msg.empty() && msg.back() == '\n')
      msg.pop_back();
    if (msg.empty())
      return;
    std::ostream &out = (level >= LogLevel::Warn) ? std::cerr : std::cout;
    out << color(level) << prefix(level) << msg << defaultColor << '\n';
  }

  template <typename T> LogStream &operator<<(T &&val) {
    if (active)
      oss << std::forward<T>(val);
    return *this;
  }
};

class Logger {
public:
  static LogLevel minLevel;
  static bool enabled;

  static LogStream debug() {
    return {LogLevel::Debug, enabled && LogLevel::Debug >= minLevel};
  }
  static LogStream info() {
    return {LogLevel::Info, enabled && LogLevel::Info >= minLevel};
  }
  static LogStream warn() {
    return {LogLevel::Warn, enabled && LogLevel::Warn >= minLevel};
  }
  static LogStream error() {
    return {LogLevel::Error, enabled && LogLevel::Error >= minLevel};
  }
};

inline LogLevel Logger::minLevel = LogLevel::Debug;
inline bool Logger::enabled = true;

#endif

} // namespace oILAB

#endif

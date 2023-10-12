bool compare_string_vectors(const std::vector<std::string>& vector1, const std::vector<std::string>& vector2) {
  if (vector1.size() != vector2.size()) {
    return false;
  }
  for (size_t i = 0; i < vector1.size(); ++i) {
    if (vector1[i] != vector2[i]) {
      return false;
    }
  }
  return true;
}







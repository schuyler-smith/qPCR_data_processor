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





void qPCR_input_headers(const std::vector<std::string>& input_files) {
  std::vector<std::string> headers;
  for (std::string file : input_files) {
    std::vector<std::string> file_headers;
    file_headers = read_csv_headers(file);
    if (!contained_in_vector(headers, {ct_col, negative_control_col, standard_id_col, non_template_control})) {
      std::cerr << "Input file " << file << " missing required input fields." << std::endl;
    }
  }
}

std::vector<std::string> read_csv_headers(const std::string& file_name) {
  std::string csv_file_name = file_name;
  std::ifstream csv_file(csv_file_name);
  if (!csv_file.is_open()) {
    std::cerr << "Failed to open the CSV file." << std::endl;
  }
  std::string header_line;
  std::vector<std::string> columns;
  std::string column;
  if (std::getline(csv_file, header_line)) {
    std::istringstream header_stream(header_line);
    while (std::getline(header_stream, column, ',')) {
      columns.push_back(column);
    }
  } else {
    std::cerr << "Failed to read the header line from the CSV file." << std::endl;
  }
  csv_file.close();
  return(columns);
}

bool contained_in_vector(const std::vector<std::string>& vector1, const std::vector<std::string>& vector2) {
  for (const std::string& arg : vector2) {
    if (std::find(vector1.begin(), vector1.end(), arg) == vector1.end()) {
      return false;
    }
  }
  return true;
}
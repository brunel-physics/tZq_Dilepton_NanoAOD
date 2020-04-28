


void findif(){

  std::vector<int> v{1, 1, 0, 0, 0, 0, 1, 0, 1};

  std::vector<size_t> results;
  floats check{};

  auto it = std::find_if(std::begin(v), std::end(v), [](int i){return i == 1;});
  while (it != std::end(v)) {
     results.emplace_back(std::distance(std::begin(v), it));
     it = std::find_if(std::next(it), std::end(v), [](int i){return i == 1;});
  }

  for(int i = 0; i < results.size(); i++){
  	check.push_back( tightjets.at( results.at(i) ) );
  }

  return check;

}

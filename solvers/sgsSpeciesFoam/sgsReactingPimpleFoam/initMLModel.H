// model for the concentration face value
std::string c_bulk_model_name = "mlModels/c_bulk_model.pt";
std::shared_ptr<torch::jit::script::Module> c_bulk_module = torch::jit::load(c_bulk_model_name);
assert(c_bulk_module != nullptr);
std::cout << "Module for face value ok.\n";

// model for the concentration gradient at the interface
std::string dc_sig_model_name = "mlModels/dc_sig_model.pt";
std::shared_ptr<torch::jit::script::Module> dc_sig_module = torch::jit::load(dc_sig_model_name);
assert(dc_sig_module != nullptr);
std::cout << "Module for interface gradient ok.\n";

// verification
torch::Tensor testTensor = torch::ones({1, 3});
testTensor[0][0] = 0.9515819615;
testTensor[0][1] = 0.2237687949;
testTensor[0][2] = 0.5;
std::vector<torch::jit::IValue> inputTest{testTensor};
std::cout << "Input tensor: " << inputTest << std::endl;
torch::Tensor outputTest = dc_sig_module->forward(inputTest).toTensor();
std::cout << "Result tensor: " << outputTest << std::endl;

// model for the concentration gradient at the interface
std::string dc_bulk_model_name = "mlModels/dc_bulk_model.pt";
std::shared_ptr<torch::jit::script::Module> dc_bulk_module = torch::jit::load(dc_bulk_model_name);
assert(dc_bulk_module != nullptr);
std::cout << "Module for bulk gradient ok.\n";

// feature tensor with dimension [nInterfaceCells x nInputFeatures]
int nInputFeatures = 3;
torch::Tensor featureTensor = torch::ones({nInterfaceCells, nInputFeatures});

IOdictionary mlModelProperties
(
    IOobject
    (
        "mlModelProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

std::vector<double> featureMinValues{};
std::vector<double> featureMaxValues{};
std::vector<double> labelMinValues{};
std::vector<double> labelMaxValues{};

// feature 1 - cav (average concentration value)
featureMinValues.push_back(dimensionedScalar(mlModelProperties.lookup("cavMin")).value());
featureMaxValues.push_back(dimensionedScalar(mlModelProperties.lookup("cavMax")).value());

// feature 2 - rf (distance of first face from interface)
featureMinValues.push_back(dimensionedScalar(mlModelProperties.lookup("rfMin")).value());
featureMaxValues.push_back(dimensionedScalar(mlModelProperties.lookup("rfMax")).value());

// feature 3 - dam (Damkoehler number)
featureMinValues.push_back(dimensionedScalar(mlModelProperties.lookup("damMin")).value());
featureMaxValues.push_back(dimensionedScalar(mlModelProperties.lookup("damMax")).value());

// label 1 - cbulk (concentration values at first face)
labelMinValues.push_back(dimensionedScalar(mlModelProperties.lookup("cbulkMin")).value());
labelMaxValues.push_back(dimensionedScalar(mlModelProperties.lookup("cbulkMax")).value());

// label 2 - dcbulk (concentration gradient at first face)
labelMinValues.push_back(dimensionedScalar(mlModelProperties.lookup("dcbulkMin")).value());
labelMaxValues.push_back(dimensionedScalar(mlModelProperties.lookup("dcbulkMax")).value());

// label 3 - dcsig (concentration gradient at interface)
labelMinValues.push_back(dimensionedScalar(mlModelProperties.lookup("dcsigMin")).value());
labelMaxValues.push_back(dimensionedScalar(mlModelProperties.lookup("dcsigMax")).value());

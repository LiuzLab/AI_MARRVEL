// modules/model_interpreter/main.nf
process modelInterpreterShapPy {
    cpus '48'
    publishDir "./published${context.publishDirSuffix}"
    
    input:
    path modelPath
    path inputDataPath
    val context
    
    output:
    path 'shapOutputs/shap_values.json', emit: shapValuesPath
    
    """
    mkdir -p ./shapOutputs
    python ${moduleDir}/main.py \
        --model-path ${modelPath} \
        --input-path ${inputDataPath} \
        --output-path ./shapOutputs/shap_values.json
    """
}

workflow modelInterpreter {
    take:
        modelPath
        inputDataPath
        context
    main:
        modelInterpreterShapPy(modelPath, inputDataPath, context)
    emit:
        shapValuesPath = modelInterpreterShapPy.out.shapValuesPath
}

workflow {
    // modelInterpreter()
}
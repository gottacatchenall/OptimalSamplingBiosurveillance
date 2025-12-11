using SpeciesDistributionToolkit
using DataFrames, CSV
using Statistics
using EvoTrees
using JSON
using CairoMakie
const SDT = SpeciesDistributionToolkit


function group_occurrences_by_species(occurrences, minimum_occurrences=50)
    function extract_binomial_name(species_string)
        name_parts = split(species_string, " ")[1:2]
        return name_parts[1] * " " * name_parts[2]
    end
    
    unique_species = unique(map(x -> extract_binomial_name(x.what), occurrences))
    
    species_dict = Dict()
    for species_name in unique_species
        species_occurrences = Occurrences(
            filter(x -> extract_binomial_name(x.what) == species_name, SDT.elements(occurrences))
        )
        if length(species_occurrences) > minimum_occurrences
            species_dict[species_name] = species_occurrences
        end
    end
    
    return species_dict
end

function generate_pseudoabsences(
    presence_layer, 
    buffer_distance_km, 
    class_balance_ratio
)
    background = pseudoabsencemask(DistanceToEvent, presence_layer)
    buffer = pseudoabsencemask(WithinRadius, presence_layer; distance = buffer_distance_km)
    
    nodata_mask = nodata(buffer, true) 
    background.indices .= nodata_mask.indices
    
    num_pseudoabsences = Int(round(class_balance_ratio * sum(presence_layer)))
    return backgroundpoints(background, num_pseudoabsences)
end

function prepare_training_data(layers, presence_layer, absence_layer)
    X = Matrix(hcat([
        vcat(layer[findall(presence_layer)], layer[findall(absence_layer)]) 
        for layer in layers
    ]...)')
    
    y = Bool.(vcat(
        [1 for _ in findall(presence_layer)], 
        [0 for _ in findall(absence_layer)]
    ))
    
    return X, y
end

function train_model(X_train, y_train)
    return EvoTrees.fit(
        EvoTreeGaussian(),
        x_train = X_train',
        y_train = y_train,
    )
end

function predict_distribution(model, feature_matrix)
    return EvoTrees.predict(model, feature_matrix')
end

function create_prediction_layer(model, environmental_layers)
    prediction_layer = deepcopy(environmental_layers[begin])
    uncertainty_layer = deepcopy(environmental_layers[begin])
    
    feature_matrix = Matrix(hcat([
        [layer[i] for layer in environmental_layers] 
        for i in eachindex(environmental_layers[1])
    ]...))
    
    predictions = predict_distribution(model, feature_matrix)
    
    prediction_layer.grid[findall(prediction_layer.indices)] .= predictions[:, 1]
    uncertainty_layer.grid[findall(prediction_layer.indices)] .= predictions[:, 2]
    
    return prediction_layer, uncertainty_layer
end

function calculate_evaluation_metrics(y_true, y_predicted, thresholds=0:0.001:1)
    # Calculate confusion matrices across all thresholds
    confusion_matrices = [ConfusionMatrix(y_predicted .> t, y_true) for t in thresholds]
    false_positive_rates, true_positive_rates = fpr.(confusion_matrices), tpr.(confusion_matrices)
    
    # Calculate ROC-AUC using trapezoidal rule
    roc_dx = [reverse(false_positive_rates)[i] - reverse(false_positive_rates)[i-1] for i in 2:length(false_positive_rates)]
    roc_dy = [reverse(true_positive_rates)[i] + reverse(true_positive_rates)[i-1] for i in 2:length(true_positive_rates)]
    roc_auc = sum(roc_dx .* (roc_dy ./ 2.0))
    
    # Calculate PR-AUC using trapezoidal rule
    precisions = ppv.(confusion_matrices)
    pr_dx = [reverse(true_positive_rates)[i] - reverse(true_positive_rates)[i-1] for i in 2:length(true_positive_rates)]
    pr_dy = [reverse(precisions)[i] + reverse(precisions)[i-1] for i in 2:length(precisions)]
    pr_auc = sum(pr_dx .* (pr_dy ./ 2.0))
    
    # Find optimal threshold using True Skill Statistic
    _, threshold_index = findmax(trueskill.(confusion_matrices))
    optimal_threshold = thresholds[threshold_index]
    
    return Dict(
        :prauc => pr_auc,
        :rocauc => roc_auc,
        :tss => trueskill(confusion_matrices[threshold_index]),
        :mcc => mcc(confusion_matrices[threshold_index]),
        :threshold => optimal_threshold
    )
end

function aggregate_fold_statistics(fold_stats)
    return Dict(
        metric => Dict("mean" => mean(values), "std" => std(values))
        for metric in keys(first(fold_stats)) 
        for values in [[fold[metric] for fold in fold_stats]]
    )
end

function fit_sdm(
    occurrences,
    environmental_layers;
    pseudoabsence_buffer_distance = 25.0,
    class_balance = 1.0,
    k = 4
)
    presence_layer = mask(environmental_layers[begin], occurrences)
    absence_layer = generate_pseudoabsences(presence_layer, pseudoabsence_buffer_distance, class_balance)
    
    features, labels = prepare_training_data(environmental_layers, presence_layer, absence_layer)
    fold_indices = SDeMo.kfold(labels, features, k=k)
    
        
    # Storage for results
    true_labels = Bool[]
    out_of_fold_predictions = Float32[]

    # Train and evaluate each fold
    for (train_idx, validation_idx) in fold_indices
        model = train_model(features[:, train_idx], labels[train_idx])
        
        # Evaluate on validation set
        validation_predictions = predict_distribution(model, features[:, validation_idx])[:, 1]

        true_labels = vcat(true_labels, labels[validation_idx])
        out_of_fold_predictions = vcat(out_of_fold_predictions, validation_predictions)
    end
    
    fit_stats = calculate_evaluation_metrics(true_labels, out_of_fold_predictions)

    model = train_model(features, labels)
    prediction, uncertainty = create_prediction_layer(model, environmental_layers)    
    
    return Dict(
        :prediction => prediction, 
        :uncertainty => uncertainty, 
        :metrics => fit_stats, 
        :presences => presence_layer, 
        :absences => absence_layer
    )
end


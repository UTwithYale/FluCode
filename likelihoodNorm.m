function likelihood = likelihoodNorm(sigma, data)
    likelihood = -sum(log(normpdf(data, 0, sigma)));
    if isinf(likelihood)
        likelihood = 10^10;
    end
end
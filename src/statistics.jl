mutable struct Averager
    sum :: Float64
    count :: Int64

    function Averager()
        return new(0.0, 0)
    end
end

function add_sample!(acc::Averager, sample)
    acc.sum += sample
    acc.count += 1
end

mean(acc::Averager) = acc.sum / acc.count

function clear!(acc::Averager)
    acc.sum = 0.0
    acc.count = 0
end


mutable struct Covariance
    sum_x  :: Float64
    sum_y  :: Float64
    sum_xy :: Float64
    count  :: Int64

    function Covariance()
        return new(0.0, 0.0, 0.0, 0)
    end
end

function add_sample!(acc::Covariance, x, y)
    acc.sum_x  += x
    acc.sum_y  += y
    acc.sum_xy += x * y
    acc.count  += 1
end

function covariance(acc::Covariance)
    n = acc.count
    n < 2 && return 0.0  # Covariance undefined for n < 2
    
    mean_x  = acc.sum_x / n
    mean_y  = acc.sum_y / n
    mean_xy = acc.sum_xy / n
    
    # Sample covariance (divides by n - 1)
    return (mean_xy - mean_x * mean_y) * n / (n - 1)
end

function clear!(acc::Covariance)
    acc.sum = 0.0
    acc.count = 0
end
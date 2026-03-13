mutable struct Averager
    sum :: Float64
    count :: Int32

    function Averager()
        return new(0.0, 0)
    end
end

function add_sample!(accumulator :: Averager, sample)
    accumulator.sum += sample
    accumulator.count += 1
end

mean(accumulator :: Averager) = accumulator.sum / accumulator.count
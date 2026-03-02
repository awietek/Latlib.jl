
# have a function that checks if a number is close to an integer within a certain tolerance
is_whole(x; atol=1e-8) = abs(x - round(x)) ≤ atol
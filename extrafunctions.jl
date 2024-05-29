function check_values_in_range(vector::Vector)
    if all(0 <= x <= 1 for x in vector)
        #println("Todos los valores están en el rango [0,1].")
        return true
    else
        #println("No todos los valores están en el rango [0,1].")
        return false
    end
end
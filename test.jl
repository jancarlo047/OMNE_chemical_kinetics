function check_values_in_range(vector)
    if all(0 <= x <= 1 for x in vector)
        println("Todos los valores están en el rango [0,1].")
    else
        println("No todos los valores están en el rango [0,1].")
    end
end

# Prueba la función con tu vector
vector = [0.1, 0.2, 0.3, 0.4, 0.5]
check_values_in_range(vector)
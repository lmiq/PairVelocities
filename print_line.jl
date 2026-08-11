using Printf


function colorize_number(x)
    s = @sprintf("%12.5f", x)
    color = x < 0 ? "\e[32m" : "\e[31m"
    return "$(color)$s\e[0m"
end

print_line(n, time, prev, previous_version) =
    println(@sprintf("%10i | current: %12.5f | %10s: %12.5f | %s / %s%1s",
        n, time, previous_version, prev,
        colorize_number(time - prev),
        colorize_number(100*(time-prev)/prev),
        "%"
    ))


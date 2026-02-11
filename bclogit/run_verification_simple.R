# Open log file
log_con <- file("verification_log.txt", open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

library(devtools)

tryCatch(
    {
        load_all(".", quiet = FALSE)
    },
    error = function(e) {
        print(e)
        quit(status = 1)
    }
)

# Mock data
set.seed(123)
n <- 20
data <- data.frame(
    y = rbinom(n, 1, 0.5),
    x1 = rnorm(n),
    x2 = rnorm(n),
    trt = rbinom(n, 1, 0.5),
    id = rep(1:(n / 2), each = 2)
)

print("Attempting bclogit without treatment...")
tryCatch(
    {
        bclogit(y ~ x1 + x2, data = data, strata = id)
        print("FAILURE: bclogit DID NOT throw an error!")
    },
    error = function(e) {
        print(paste("SUCCESS: bclogit threw error:", e$message))
    }
)

print("Attempting bclogit.default without treatment...")
X <- model.matrix(~ x1 + x2, data = data)
tryCatch(
    {
        bclogit.default(response = data$y, data = X, strata = data$id)
        print("FAILURE: bclogit.default DID NOT throw an error!")
    },
    error = function(e) {
        print(paste("SUCCESS: bclogit.default threw error:", e$message))
    }
)

# Close sinks
sink(type = "message")
sink(type = "output")
close(log_con)

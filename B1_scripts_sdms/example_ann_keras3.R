install.packages("keras3")
# load("keras3")
# library("keras3")
# remotes::install_github("rstudio/keras")
reticulate::install_python(version = '<version>')
# reticulate::install_python()
library(reticulate)
# virtualenv_dir <- "C:/Users/cholleem/.virtualenvs"
# dir.create(virtualenv_dir, showWarnings = FALSE)
# use_virtualenv(virtualenv_dir, required = TRUE)

Sys.getenv()
# Sys.setenv(HOME = "C:/Users/cholleem")
# Sys.setenv(R_USER = "C:/Users/cholleem")
# 
# "C:\Users\cholleem"


py_config()
# virtualenv_create("r-keras", python = "C:/Users/cholleem/AppData/Local/r-reticulate/r-reticulate/pyenv/pyenv-win/versions/3.10.11/python.exe")
# 
keras3::install_keras(backend = "tensorflow")
keras3::install_keras()

library(keras3)
mnist <- dataset_mnist()
x_train <- mnist$train$x
y_train <- mnist$train$y
x_test <- mnist$test$x
y_test <- mnist$test$y

# reshape
x_train <- array_reshape(x_train, c(nrow(x_train), 784))
x_test <- array_reshape(x_test, c(nrow(x_test), 784))
# rescale
x_train <- x_train / 255
x_test <- x_test / 255

y_train <- to_categorical(y_train, 10)
y_test <- to_categorical(y_test, 10)

model <- keras_model_sequential(input_shape = c(784))
model |>
    layer_dense(units = 256, activation = 'relu') |>
    layer_dropout(rate = 0.4) |>
    layer_dense(units = 128, activation = 'relu') |>
    layer_dropout(rate = 0.3) |>
    layer_dense(units = 10, activation = 'softmax')

summary(model)

# plot(model)

model |> compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
)

history <- model |> fit(
    x_train, y_train,
    epochs = 30, batch_size = 128,
    validation_split = 0.2
)

plot(history)

model |> evaluate(x_test, y_test)

probs <- model |> predict(x_test)

model |> save_model('example_model.keras')

model <- load_model('example_model.keras')

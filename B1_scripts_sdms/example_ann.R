model <- keras_model_sequential()

model %>%
    
    # Adds a densely-connected layer with 64 units to the model:
    layer_dense(units = 64, activation = 'relu') %>%
    
    # Add another:
    layer_dense(units = 64, activation = 'relu') %>%
    
    # Add a softmax layer with 10 output units:
    layer_dense(units = 10, activation = 'softmax')

model %>% compile(
    optimizer = 'adam',
    loss = 'categorical_crossentropy',
    metrics = list('accuracy')
)

data <- matrix(rnorm(1000 * 32), nrow = 1000, ncol = 32)
labels <- matrix(rnorm(1000 * 10), nrow = 1000, ncol = 10)

model %>% fit(
    data,
    labels,
    epochs = 10,
    batch_size = 32
)


# model %>% save_model_weights_tf('my_model/')
# Save entire model to the SavedModel format
model %>% save_model_tf('my_model/')

# Recreate the exact same model, including weights and optimizer.
model1 <- load_model_tf('my_model/')

pred1 <- model1 %>% predict(data, batch_size = 32)
pred <- model %>% predict(data, batch_size = 32)

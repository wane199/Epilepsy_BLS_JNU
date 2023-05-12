## Deep Learning and Scientific Computing with R torch(https://skeydan.github.io/Deep-Learning-and-Scientific-Computing-with-R-torch/)
# 19  Image segmentation
# 19.1 Segmentation vs. classification
# 19.2 U-Net, a “classic” in image segmentation
# 19.3 U-Net – a torch implementation
# 19.3.1 Encoder
rm(list = ls())
library(torch)
library(torchvision)
library(luz)

model <- model_mobilenet_v2(pretrained = TRUE)
model

model$features

model$features[2]$`0`$conv

encoder <- nn_module(
  initialize = function() {
    model <- model_mobilenet_v2(pretrained = TRUE)
    self$stages <- nn_module_list(list(
      nn_identity(),
      model$features[1:2],
      model$features[3:4],
      model$features[5:7],
      model$features[8:14],
      model$features[15:18]
    ))
    for (par in self$parameters) {
      par$requires_grad_(FALSE)
    }
  },
  forward = function(x) {
    features <- list()
    for (i in 1:length(self$stages)) {
      x <- self$stages[[i]](x)
      features[[length(features) + 1]] <- x
    }
    features
  }
)

sample <- torch_randn(1, 3, 224, 224)
sample_features <- encoder()(sample)
purrr::map(sample_features, purrr::compose(dim, as.array))

decoder_block <- nn_module(
  initialize = function(in_channels,
                        skip_channels,
                        out_channels) {
    self$upsample <- nn_conv_transpose2d(
      in_channels = in_channels,
      out_channels = out_channels,
      kernel_size = 2,
      stride = 2
    )
    self$activation <- nn_relu()
    self$conv <- nn_conv2d(
      in_channels = out_channels + skip_channels,
      out_channels = out_channels,
      kernel_size = 3,
      padding = "same"
    )
  },
  forward = function(x, skip) {
    x <- x %>%
      self$upsample() %>%
      self$activation()
    input <- torch_cat(list(x, skip), dim = 2)
    input %>%
      self$conv() %>%
      self$activation()
  }
)

img <- torch_randn(1, 1, 5, 5)

conv <- nn_conv2d(
  in_channels = 1,
  out_channels = 1,
  kernel_size = 3,
  stride = 2,
  padding = 1
)

convolved <- conv(img)
convolved

transposed_conv <- nn_conv_transpose2d(
  in_channels = 1,
  out_channels = 1,
  kernel_size = 3,
  stride = 2,
  padding = 1
)

upsampled <- transposed_conv(convolved)
upsampled

first_decoder_block <- decoder_block(
  in_channels = 320,
  skip_channels = 96,
  out_channels = 256
)

first_decoder_block(
  sample_features[[6]],
  sample_features[[5]]
) %>%
  dim()

decoder <- nn_module(
  initialize = function(
    decoder_channels = c(256, 128, 64, 32, 16),
    encoder_channels = c(16, 24, 32, 96, 320)) {
    encoder_channels <- rev(encoder_channels)
    skip_channels <- c(encoder_channels[-1], 3)
    in_channels <- c(encoder_channels[1], decoder_channels)
    
    depth <- length(encoder_channels)
    
    self$blocks <- nn_module_list()
    for (i in seq_len(depth)) {
      self$blocks$append(decoder_block(
        in_channels = in_channels[i],
        skip_channels = skip_channels[i],
        out_channels = decoder_channels[i]
      ))
    }
  },
  forward = function(features) {
    features <- rev(features)
    x <- features[[1]]
    for (i in seq_along(self$blocks)) {
      x <- self$blocks[[i]](x, features[[i + 1]])
    }
    x
  }
)

model <- nn_module(
  initialize = function() {
    self$encoder <- encoder()
    self$decoder <- decoder()
    self$output <- nn_conv2d(
      in_channels = 16,
      out_channels = 3,
      kernel_size = 3,
      padding = "same"
    )
  },
  forward = function(x) {
    x %>%
      self$encoder() %>%
      self$decoder() %>%
      self$output()
  }
)


# 19.4 Dogs and cats
library(torchvision)
library(torchdatasets)

dir <- "~/.torch-datasets"

ds <- oxford_pet_dataset(root = dir, download = TRUE)

pet_dataset <- torch::dataset(
  inherit = oxford_pet_dataset,
  initialize = function(...,
                        size,
                        normalize = TRUE,
                        augmentation = NULL) {
    self$augmentation <- augmentation
    input_transform <- function(x) {
      x <- x %>%
        transform_to_tensor() %>%
        transform_resize(size)
      if (normalize) {
        x <- x %>%
          transform_normalize(
            mean = c(0.485, 0.456, 0.406),
            std = c(0.229, 0.224, 0.225)
          )
      }
      x
    }
    target_transform <- function(x) {
      x <- torch_tensor(x, dtype = torch_long())
      x <- x[newaxis, ..]
      # interpolation = 0 makes sure we
      # still end up with integer classes
      x <- transform_resize(x, size, interpolation = 0)
      x[1, ..]
    }
    super$initialize(
      ...,
      transform = input_transform,
      target_transform = target_transform
    )
  },
  .getitem = function(i) {
    item <- super$.getitem(i)
    if (!is.null(self$augmentation)) {
      self$augmentation(item)
    } else {
      list(x = item$x, y = item$y)
    }
  }
)




















################################
# R原生支持pytorch了，应该叫torch for R
library(torch)

# creates example tensors. x requires_grad = TRUE tells that 
# we are going to take derivatives over it.
dense <- nn_module(
  clasname = "dense",
  # the initialize function tuns whenever we instantiate the model
  initialize = function(in_features, out_features) {
    
    # just for you to see when this function is called
    cat("Calling initialize!") 
    
    # we use nn_parameter to indicate that those tensors are special
    # and should be treated as parameters by `nn_module`.
    self$w <- nn_parameter(torch_randn(in_features, out_features))
    self$b <- nn_parameter(torch_zeros(out_features))
    
  },
  # this function is called whenever we call our model on input.
  forward = function(x) {
    cat("Calling forward!")
    torch_mm(x, self$w) + self$b
  }
)

model <- dense(3, 1)

# you can get all parameters 
model$parameters

# create an input tensor
x <- torch_randn(10, 3)
y_pred <- model(x)
y_pred


########################################
# Brain image segmentation with torch(https://blogs.rstudio.com/ai/posts/2020-11-30-torch-brain-segmentation/)
# Data
# deep learning (incl. dependencies)
rm(list = ls())
library(torch)
library(torchvision)

# data wrangling
library(tidyverse)
library(zeallot)

# image processing and visualization
library(magick)
library(cowplot)

# dataset loading 
library(pins)
library(zip)

torch_manual_seed(123)
set.seed(123)

# use your own kaggle.json here
pins::board_register_kaggle(token = "/home/wane/Downloads/kaggle.json")

# files <- pins::pin_get("mateuszbuda/lgg-mri-segmentation", board = "kaggle",  extract = FALSE)

train_dir <- "/home/wane/Downloads/Brain_MRI_segmentation/data/mri_train"
valid_dir <- "/home/wane/Downloads/Brain_MRI_segmentation/data/mri_valid"

# if(dir.exists(train_dir)) unlink(train_dir, recursive = TRUE, force = TRUE)
# if(dir.exists(valid_dir)) unlink(valid_dir, recursive = TRUE, force = TRUE)
# 
# zip::unzip(files, exdir = "/home/wane/Downloads/Brain_MRI_segmentation/data")
# 
# file.rename("/home/wane/Downloads/Brain_MRI_segmentation/data/kaggle_3m", train_dir)
# 
# # this is a duplicate, again containing kaggle_3m (evidently a packaging error on Kaggle)
# # we just remove it
# unlink("/home/wane/Downloads/Brain_MRI_segmentation/data/lgg-mri-segmentation", recursive = TRUE)
# 
# dir.create(valid_dir)
# 
# patients <- list.dirs(train_dir, recursive = FALSE)
# 
# valid_indices <- sample(1:length(patients), 30)
# 
# for (i in valid_indices) {
#   dir.create(file.path(valid_dir, basename(patients[i])))
#   for (f in list.files(patients[i])) {    
#     file.rename(file.path(train_dir, basename(patients[i]), f), file.path(valid_dir, basename(patients[i]), f))    
#   }
#   unlink(file.path(train_dir, basename(patients[i])), recursive = TRUE)
# }

# Dataset
brainseg_dataset <- dataset(
  name = "brainseg_dataset",
  
  initialize = function(img_dir,
                        augmentation_params = NULL,
                        random_sampling = FALSE) {
    self$images <- tibble(
      img = grep(
        list.files(
          img_dir,
          full.names = TRUE,
          pattern = "tif",
          recursive = TRUE
        ),
        pattern = 'mask',
        invert = TRUE,
        value = TRUE
      ),
      mask = grep(
        list.files(
          img_dir,
          full.names = TRUE,
          pattern = "tif",
          recursive = TRUE
        ),
        pattern = 'mask',
        value = TRUE
      )
    )
    self$slice_weights <- self$calc_slice_weights(self$images$mask)
    self$augmentation_params <- augmentation_params
    self$random_sampling <- random_sampling
  },
  
  .getitem = function(i) {
    index <-
      if (self$random_sampling == TRUE)
        sample(1:self$.length(), 1, prob = self$slice_weights)
    else
      i
    
    img <- self$images$img[index] %>%
      image_read() %>%
      transform_to_tensor() 
    mask <- self$images$mask[index] %>%
      image_read() %>%
      transform_to_tensor() %>%
      transform_rgb_to_grayscale() %>%
      torch_unsqueeze(1)
    
    img <- self$min_max_scale(img)
    
    if (!is.null(self$augmentation_params)) {
      scale_param <- self$augmentation_params[1]
      c(img, mask) %<-% self$resize(img, mask, scale_param)
      
      rot_param <- self$augmentation_params[2]
      c(img, mask) %<-% self$rotate(img, mask, rot_param)
      
      flip_param <- self$augmentation_params[3]
      c(img, mask) %<-% self$flip(img, mask, flip_param)
      
    }
    list(img = img, mask = mask)
  },
  
  .length = function() {
    nrow(self$images)
  },
  
  calc_slice_weights = function(masks) {
    weights <- map_dbl(masks, function(m) {
      img <-
        as.integer(magick::image_data(image_read(m), channels = "gray"))
      sum(img / 255)
    })
    
    sum_weights <- sum(weights)
    num_weights <- length(weights)
    
    weights <- weights %>% map_dbl(function(w) {
      w <- (w + sum_weights * 0.1 / num_weights) / (sum_weights * 1.1)
    })
    weights
  },
  
  min_max_scale = function(x) {
    min = x$min()$item()
    max = x$max()$item()
    x$clamp_(min = min, max = max)
    x$add_(-min)$div_(max - min + 1e-5)
    x
  },
  
  resize = function(img, mask, scale_param) {
    img_size <- dim(img)[2]
    rnd_scale <- runif(1, 1 - scale_param, 1 + scale_param)
    img <- transform_resize(img, size = rnd_scale * img_size)
    mask <- transform_resize(mask, size = rnd_scale * img_size)
    diff <- dim(img)[2] - img_size
    if (diff > 0) {
      top <- ceiling(diff / 2)
      left <- ceiling(diff / 2)
      img <- transform_crop(img, top, left, img_size, img_size)
      mask <- transform_crop(mask, top, left, img_size, img_size)
    } else {
      img <- transform_pad(img,
                           padding = -c(
                             ceiling(diff / 2),
                             floor(diff / 2),
                             ceiling(diff / 2),
                             floor(diff / 2)
                           ))
      mask <- transform_pad(mask, padding = -c(
        ceiling(diff / 2),
        floor(diff /
                2),
        ceiling(diff /
                  2),
        floor(diff /
                2)
      ))
    }
    list(img, mask)
  },
  
  rotate = function(img, mask, rot_param) {
    rnd_rot <- runif(1, 1 - rot_param, 1 + rot_param)
    img <- transform_rotate(img, angle = rnd_rot)
    mask <- transform_rotate(mask, angle = rnd_rot)
    
    list(img, mask)
  },
  
  flip = function(img, mask, flip_param) {
    rnd_flip <- runif(1)
    if (rnd_flip > flip_param) {
      img <- transform_hflip(img)
      mask <- transform_hflip(mask)
    }
    
    list(img, mask)
  }
)

train_ds <- brainseg_dataset(
  train_dir,
  augmentation_params = c(0.05, 15, 0.5),
  random_sampling = TRUE
)

length(train_ds)

valid_ds <- brainseg_dataset(
  valid_dir,
  augmentation_params = NULL,
  random_sampling = FALSE
)

length(valid_ds)

par(mfrow = c(1, 2), mar = c(0, 1, 0, 1))

img_and_mask <- valid_ds[27]
img <- img_and_mask[[1]]
mask <- img_and_mask[[2]]

img$permute(c(2, 3, 1)) %>% as.array() %>% as.raster() %>% plot()
mask$squeeze() %>% as.array() %>% as.raster() %>% plot()

img_and_mask <- valid_ds[77]
img <- img_and_mask[[1]]
mask <- img_and_mask[[2]]

imgs <- map (1:24, function(i) {
  
  # scale factor; train_ds really uses 0.05
  c(img, mask) %<-% valid_ds$resize(img, mask, 0.2) 
  c(img, mask) %<-% valid_ds$flip(img, mask, 0.5)
  # rotation angle; train_ds really uses 15
  c(img, mask) %<-% valid_ds$rotate(img, mask, 90) 
  img %>%
    transform_rgb_to_grayscale() %>%
    as.array() %>%
    as_tibble() %>%
    rowid_to_column(var = "Y") %>%
    gather(key = "X", value = "value", -Y) %>%
    mutate(X = as.numeric(gsub("V", "", X))) %>%
    ggplot(aes(X, Y, fill = value)) +
    geom_raster() +
    theme_void() +
    theme(legend.position = "none") +
    theme(aspect.ratio = 1)
  
})

plot_grid(plotlist = imgs, nrow = 4)

batch_size <- 4
train_dl <- dataloader(train_ds, batch_size)
valid_dl <- dataloader(valid_ds, batch_size)

# Model
unet <- nn_module(
  "unet",
  
  initialize = function(channels_in = 3,
                        n_classes = 1,
                        depth = 5,
                        n_filters = 6) {
    
    self$down_path <- nn_module_list()
    
    prev_channels <- channels_in
    for (i in 1:depth) {
      self$down_path$append(down_block(prev_channels, 2 ^ (n_filters + i - 1)))
      prev_channels <- 2 ^ (n_filters + i -1)
    }
    
    self$up_path <- nn_module_list()
    
    for (i in ((depth - 1):1)) {
      self$up_path$append(up_block(prev_channels, 2 ^ (n_filters + i - 1)))
      prev_channels <- 2 ^ (n_filters + i - 1)
    }
    
    self$last = nn_conv2d(prev_channels, n_classes, kernel_size = 1)
  },
  
  forward = function(x) {
    
    blocks <- list()
    
    for (i in 1:length(self$down_path)) {
      x <- self$down_path[[i]](x)
      if (i != length(self$down_path)) {
        blocks <- c(blocks, x)
        x <- nnf_max_pool2d(x, 2)
      }
    }
    
    for (i in 1:length(self$up_path)) {  
      x <- self$up_path[[i]](x, blocks[[length(blocks) - i + 1]]$to(device = device))
    }
    
    torch_sigmoid(self$last(x))
  }
)

down_block <- nn_module(
  "down_block",
  
  initialize = function(in_size, out_size) {
    self$conv_block <- conv_block(in_size, out_size)
  },
  
  forward = function(x) {
    self$conv_block(x)
  }
)

up_block <- nn_module(
  "up_block",
  
  initialize = function(in_size, out_size) {
    
    self$up = nn_conv_transpose2d(in_size,
                                  out_size,
                                  kernel_size = 2,
                                  stride = 2)
    self$conv_block = conv_block(in_size, out_size)
  },
  
  forward = function(x, bridge) {
    
    up <- self$up(x)
    torch_cat(list(up, bridge), 2) %>%
      self$conv_block()
  }
)

conv_block <- nn_module( 
  "conv_block",
  
  initialize = function(in_size, out_size) {
    
    self$conv_block <- nn_sequential(
      nn_conv2d(in_size, out_size, kernel_size = 3, padding = 1),
      nn_relu(),
      nn_dropout(0.6),
      nn_conv2d(out_size, out_size, kernel_size = 3, padding = 1),
      nn_relu()
    )
  },
  
  forward = function(x){
    self$conv_block(x)
  }
)

device <- torch_device(if(cuda_is_available()) "cuda" else "cpu")
model <- unet(depth = 5)$to(device = device)

# Optimization
calc_dice_loss <- function(y_pred, y_true) {
  
  smooth <- 1
  y_pred <- y_pred$view(-1)
  y_true <- y_true$view(-1)
  intersection <- (y_pred * y_true)$sum()
  
  1 - ((2 * intersection + smooth) / (y_pred$sum() + y_true$sum() + smooth))
}

dice_weight <- 0.3

optimizer <- optim_sgd(model$parameters, lr = 0.1, momentum = 0.9)

num_epochs <- 20

scheduler <- lr_one_cycle(
  optimizer,
  max_lr = 0.1,
  steps_per_epoch = length(train_dl),
  epochs = num_epochs
)

# Training
train_batch <- function(b) {
  
  optimizer$zero_grad()
  output <- model(b[[1]]$to(device = device))
  target <- b[[2]]$to(device = device)
  
  bce_loss <- nnf_binary_cross_entropy(output, target)
  dice_loss <- calc_dice_loss(output, target)
  loss <-  dice_weight * dice_loss + (1 - dice_weight) * bce_loss
  
  loss$backward()
  optimizer$step()
  scheduler$step()
  
  list(bce_loss$item(), dice_loss$item(), loss$item())
  
}

valid_batch <- function(b) {
  
  output <- model(b[[1]]$to(device = device))
  target <- b[[2]]$to(device = device)
  
  bce_loss <- nnf_binary_cross_entropy(output, target)
  dice_loss <- calc_dice_loss(output, target)
  loss <-  dice_weight * dice_loss + (1 - dice_weight) * bce_loss
  
  list(bce_loss$item(), dice_loss$item(), loss$item())
  
}

for (epoch in 1:num_epochs) {
  
  model$train()
  train_bce <- c()
  train_dice <- c()
  train_loss <- c()
  
  coro::loop(for (b in train_dl) {
    c(bce_loss, dice_loss, loss) %<-% train_batch(b)
    train_bce <- c(train_bce, bce_loss)
    train_dice <- c(train_dice, dice_loss)
    train_loss <- c(train_loss, loss)
  })
  
  torch_save(model, paste0("model_", epoch, ".pt"))
  
  cat(sprintf("\nEpoch %d, training: loss:%3f, bce: %3f, dice: %3f\n",
              epoch, mean(train_loss), mean(train_bce), mean(train_dice)))
  
  model$eval()
  valid_bce <- c()
  valid_dice <- c()
  valid_loss <- c()
  
  i <- 0
  coro::loop(for (b in tvalid_dl) {
    
    i <<- i + 1
    c(bce_loss, dice_loss, loss) %<-% valid_batch(b)
    valid_bce <- c(valid_bce, bce_loss)
    valid_dice <- c(valid_dice, dice_loss)
    valid_loss <- c(valid_loss, loss)
    
  })
  
  cat(sprintf("\nEpoch %d, validation: loss:%3f, bce: %3f, dice: %3f\n",
              epoch, mean(valid_loss), mean(valid_bce), mean(valid_dice)))
}

# Evaluation
saved_model <- torch_load("model_20.pt") 

model <- saved_model
model$eval()

saved_model <- torch_load("model_20.pt") 

model <- saved_model
model$eval()


















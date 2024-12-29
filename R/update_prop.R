#
# update_prop(ds, target, sources, set, get, combine = function(l) Reduce('+', l), override = function(ds, target, v) v)
#
#     ds:       dataset with properties to be updated
#     target:   key of target object in dataset
#     sources:  keys of source objects in dataset
#     set:      set method for property, called as set(ds, id, value)
#     get:      get method for property, called as get(ds, id)
#     combine:  combining operator, called as combine(values)
#     override: override value function, called as override(ds, target, value)
#

update_prop <- function(ds, target, sources, set, get,
                        combine = function(l) Reduce('+', l),
                        override = function(ds, target, v) v) {
  if (length(sources) > 0) {
    av <- Map(f = function(s) get(ds, s), sources)
    set(ds, target, override(ds, target, combine(av)))
  } else
    ds
}

# helper methods for data frames

# get keys from data frame

df_get_keys <- function(df, key) df[, key]

# get ids from data frame (key="id")

df_get_ids <- function(df) df_get_keys(df, "id")

# get property by key from data frame

df_get_by_key <- function(df, key, r, c) df[df[, key] == r, c]

# get property by key="id" from data frame

df_get_by_id <- function(df, r, c) df_get_by_key(df, "id", r, c)

# set property by key in data frame

df_set_by_key <- function(df, key, r, c, v) {
  df[df[, key] == r, c] <- v
  df
}

# set property by key="id" in data frame

df_set_by_id <- function(df, r, c, v) {
  df_set_by_key(df, "id", r, c, v)
}

# 
# update_df_prop_by_id(df, target, sources, prop, ...)
#
#     df:       data frame with properties to be updated
#     target:   id of target object in dataset
#     sources:  ids of source objects in dataset
#     prop:     property in dataset to be combined and updated

update_df_prop_by_id <- function(df, target, sources, prop, ...) {
  update_prop(df, target, sources,
              function(df, id, v) df_set_by_id(df, id, prop, v),
              function(df, id, v) df_get_by_id(df, id, prop),
              ...
  )
}

# validate tree, table

validate_df_with_id <- function(tree, df, prop, ...) {
  validate_table(tree, df, df_get_ids, function(df, r) df_get_by_id(df, r, prop), ...)
}

# combine operator for root-sum-square

combine_rss <- function(vl) {
  sqrt(Reduce(f = '+', x = Map(f = function(x) x^2, x = vl)))
}


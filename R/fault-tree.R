df_get_fault_props <- function(df, id) {
  list(
    type = df_get_by_id(df, id, "type"),
    prob = df_get_by_id(df, id, "prob")
  )
}

df_set_fault_props <- function(df, id, v) {
  df_set_by_id(df, id, "prob", v$prob)
}

override_fault_props <- function(df, id, vl) {
  op <- if (df_get_by_id(df, id, "type") == "and") "*" else "+"
  list(
    prob = Reduce(
      f = op,
      Map(f = function(v) v$prob, vl)
    )
  )
}

update_fault_props <- function(ds, parent_key, child_keys) {
  update_prop(
    ds,
    target = parent_key,
    sources = child_keys,
    set = df_set_fault_props,
    get = df_get_fault_props,
    combine = function(vl) vl,
    override = override_fault_props
  )
}

validate_fault_props <- function(fp) {
  if (fp$type != "basic") stop(sprintf("invalid leaf node type %s", fp$type))
  if (!is.numeric(fp$prob) || fp$prob < 0.0 || fp$prob > 1.0) stop(sprintf("invalid probability value %f", fp$prob))
  TRUE
}

validate_fault_props_table <- function(tree, df) {
  validate_table(tree, df, df_get_ids, df_get_fault_props, validate_fault_props)
}

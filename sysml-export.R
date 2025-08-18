library(readr)
suppressPackageStartupMessages(library(dplyr))
library(jsonlite)
library(uuid)
# suppressPackageStartupMessages(library(purrr))

library(rollupTree)
library(massProps)

# Utility Functions

sysml_get_ids <- function(ds) {
  Map(
    f = function(i) i[["declaredShortName"]],
    Filter(
      f = function(i) i[["@type"]] == "PartUsage",
      ds
    )
  ) |> unlist() |> unname()
}

sysml_get_part_with_id <- function(ds, id) {
  Filter(
    f = function(i) i[["@type"]] == "PartUsage" && i[["declaredShortName"]] == id,
    ds
  )[[1]]
}

sysml_get_relations_with_source <- function(ds, type, source) {
  Filter(
    f = function(i) i[["@type"]] == type && i[["source"]][[1]][["@id"]] == source,
    ds
  )
}

sysml_get_relations_with_target <- function(ds, type, target) {
  Filter(
    f = function(i) i[["@type"]] == type && i[["target"]][[1]][["@id"]] == target,
    ds
  )
}

sysml_get_attribute_usages_for_part <- function(ds, part, attribute) {
  Filter(
    f = function(i) i[["@type"]] == "AttributeUsage" && i[["declaredName"]] == attribute,
    Map(
      f = function(i) ds[[i[["target"]][[1]][["@id"]]]],
      sysml_get_relations_with_source(ds, "FeatureMembership", part[["@id"]])
    )
  )
}

sysml_get_negated_rational <- function(ds, operator_expression) {
  ft <- Map(
    f = function(i) ds[[unlist(i[["target"]])]],
    sysml_get_relations_with_source(ds, "ParameterMembership", operator_expression[["@id"]])
  )[[1]]
  lr <- Map(
    f = function(i) ds[[unlist(i[["target"]])]],
    sysml_get_relations_with_source(ds, "FeatureValue", ft[["@id"]])
  )[[1]]
  -lr[["value"]]
}

sysml_get_attribute_values <- function(ds, attribute_usage) {
  Map(
    f = function(i) {
      fv = ds[[i[["target"]][[1]][["@id"]]]]
      ifelse(fv[["@type"]] == "OperatorExpression", sysml_get_negated_rational(ds, fv), fv[["value"]])
    },
    sysml_get_relations_with_source(ds, "FeatureValue", attribute_usage[["@id"]])
  )
}

sysml_get_by_id <- function(ds, id, property) {
  sysml_get_attribute_values(
    ds,
    sysml_get_attribute_usages_for_part(
      ds,
      sysml_get_part_with_id(ds, id),
      property)[[1]]
  ) |> unlist() |> unname()
}

sysml_get_mass_props <- function(ds, id) {
  get_mass_props(ds, id, get_by_id = sysml_get_by_id)
}

sysml_get_mass_props_and_unc <- function(ds, id) {
  get_mass_props_and_unc(ds, id, get_by_id = sysml_get_by_id)
}

sysml_get_edgelist <- function(ds) {
  Reduce(
    f = rbind,
    x = Map(
      f = function(i) data.frame(
        child = ds[[i[["target"]][[1]][["@id"]]]][["declaredShortName"]],
        parent = ds[[i[["source"]][[1]][["@id"]]]][["declaredShortName"]]
      ),
      Filter(
        f = function(i) {
          i[["@type"]] == "FeatureMembership" &&
            ds[[i[["source"]][[1]][["@id"]]]][["@type"]] == "PartUsage" &&
            ds[[i[["target"]][[1]][["@id"]]]][["@type"]] == "PartUsage"
        },
        x = ds
      )
    )
  ) |> as.matrix()
}

sysml_create_literal <- function(value) {
  uuid <- UUIDgenerate()
  list(
    `@type` = NA,
    `@id` = uuid,
    `isConstant` = FALSE,
    `isDerived` = FALSE,
    `isImpliedIncluded` = FALSE,
    `isAbstract` = FALSE,
    `isComposite` = FALSE,
    `ownedRelationship` = list(),
    `aliasIds` = list(),
    `isSufficient` = FALSE,
    `value` = value,
    `isOrdered` = FALSE,
    `elementId` = uuid,
    `isEnd` = FALSE,
    `isUnique` = TRUE,
    `isVariable` = FALSE,
    `owner` = list(
      `@id` = NA # AttributeUsage id
    ),
    `isPortion` = FALSE,
    `owningRelationship` = list(
      `@id` = NA # FeatureValue id
    ),
    `isLibraryElement` = FALSE
  )
}

sysml_create_literal_rational <- function(value) {
  lr <- sysml_create_literal(value)
  lr[["@type"]] = "LiteralRational"
  lr
}

sysml_create_literal_boolean <- function(value) {
  lr <- sysml_create_literal(value)
  lr[["@type"]] = "LiteralBoolean"
  lr
}

sysml_create_literal_string <- function(value) {
  lr <- sysml_create_literal(value)
  lr[["@type"]] = "LiteralString"
  lr
}

sysml_create_feature_value <- function() {
  uuid <- UUIDgenerate()
  list(
    `@type` = "FeatureValue",
    `@id` = uuid,
    `isInitial` = FALSE,
    `isImpliedIncluded` = FALSE,
    `isImplied` = FALSE,
    `ownedRelationship` = list(),
    `aliasIds` = list(),
    `memberElement` = list(
      `@id` = NA # Literal id
    ),
    `source` = list(
      list(
        `@id` = NA # AttributeUsage id
      )
    ),
    `target` = list(
      list(
        `@id` = NA # Literal id
      )
    ),
    `isDefault` = FALSE,
    `ownedRelatedElement` = list(
      list(
        `@id` = NA # Literal id
      )
    ),
    `owningRelatedElement` = list(
      `@id` = NA # AttributeUsage id
    ),
    `elementId` = uuid,
    `visibility` = "public",
    `isLibraryElement` = FALSE
  )
}

sysml_create_operator_expression <- function() {
  uuid <- UUIDgenerate()
  list(
    `@type` = "OperatorExpression",
    `@id` = uuid,
    `isConstant` = FALSE,
    `operator` = "-",
    `isDerived` = FALSE,
    `isImpliedIncluded` = FALSE,
    `isAbstract` = FALSE,
    `isComposite` = FALSE,
    `ownedRelationship` = list(
      list(
        `@id` = NA # ParameterMembership id
      ),
      list(
        `@id` = NA # ReturnParameterMembership id
      )
    ),
    `aliasIds` = list(),
    `isSufficient` = FALSE,
    `isOrdered` = FALSE,
    `isEnd` = FALSE,
    `elementId` = uuid,
    `isUnique` = TRUE,
    `isVariable` = FALSE,
    `owner` = list(
      `@id` = NA # AttributeUsage id
    ),
    `isPortion` = FALSE,
    `owningRelationship` = list(
      `@id` = NA # FeatureValue id (for AttributeUsage)
    ),
    `isLibraryElement` = FALSE
  )
}

sysml_create_parameter_membership <- function() {
  id = UUIDgenerate()
  list(
    `@type` = "ParameterMembership",
    `@id` = id,
    `isImpliedIncluded` = FALSE,
    `isImplied` = FALSE,
    `ownedRelationship` = list(),
    `aliasIds` = list(),
    `memberElement` = list(
      `@id` = NA # Feature id
    ),
    `source` = list(
      list(
        `@id` = NA # OperatorExpression id
      )
    ),
    `target` = list(
      list(
        `@id` =  NA # Feature id
      )
    ),
    `ownedRelatedElement` = list(
      list(
        `@id` = NA # Feature id
      )
    ),
    `owningRelatedElement` = list(
      `@id` = NA # OperatorExpression id
    ),
    `elementId` = id,
    `memberName` = "x",
    `visibility` = "private",
    `isLibraryElement` = FALSE
  )
}

sysml_create_return_parameter_membership <- function() {
  id = UUIDgenerate()
  list(
    `@type` = "ReturnParameterMembership",
    `@id` = id,
    `isImpliedIncluded` = FALSE,
    `isImplied` = FALSE,
    `ownedRelationship` = list(),
    `aliasIds` = list(),
    `memberElement` = list(
      `@id` = "7a6bc8e0-677e-43bc-988a-bb36ed27607b"
    ),
    `source` = list(
      list(
        `@id` = "27071c56-ab7f-4a4a-bd03-3e2f4967b1d5"
      )
    ),
    `target` = list(
      list(
        `@id` = "7a6bc8e0-677e-43bc-988a-bb36ed27607b"
      )
    ),
    `owningRelatedElement` = list(
      `@id` = "27071c56-ab7f-4a4a-bd03-3e2f4967b1d5"
    ),
    `ownedRelatedElement` = list(
      list(
        `@id` = "7a6bc8e0-677e-43bc-988a-bb36ed27607b"
      )
    ),
    `elementId` = id,
    `memberName` = "result",
    `visibility` = "public",
    `isLibraryElement` = FALSE
  )
}

sysml_create_feature <- function() {
  id = UUIDgenerate()
  list(
    `@type` = "Feature",
    `@id` = id,
    `direction` = "out",
    `isConstant` = false,
    `isDerived` = false,
    `isImpliedIncluded` = false,
    `isAbstract` = false,
    `isComposite` = false,
    `ownedRelationship` = list(),
    `aliasIds` = list(),
    `isSufficient` = false,
    `isOrdered` = false,
    `elementId` = id,
    `isEnd` = false,
    `isUnique` = true,
    `isVariable` = false,
    `owner` = list(
      `@id` = "27071c56-ab7f-4a4a-bd03-3e2f4967b1d5"
    ),
    `isPortion` = false,
    `owningRelationship` = list(
      `@id` = "0a18f170-03ee-4ae2-a9b8-545a4e50db73"
    ),
    `isLibraryElement` = false
  )
}

sysml_set_by_id <- function(ds, id, property, value) {
  pa <- sysml_get_part_with_id(ds, id)
  if (is.null(pa)) stop(sprintf("no part for id %s", id))
  
  au <- sysml_get_attribute_usages_for_part(ds, pa, property)[[1]]
  if (is.null(au)) stop(sprintf("no property %s for part %s", property, id))
  au_id <- au[["@id"]]
  
  # delete any existing values
  
  for (fv in sysml_get_relations_with_source(ds, "FeatureValue", au_id)) {
    ds[[fv[["target"]][[1]][["@id"]]]] <- NULL
    ds[[fv[["@id"]]]] <- NULL
  }
  
  lv <- if (is.numeric(value)) {
    sysml_create_literal_rational(value)
  } else if (is.logical(value)) {
    sysml_create_literal_boolean(value)
  } else if (is.character(value)) {
    sysml_create_literal_string(value)
  }
  lv_id <- lv[["@id"]]
  
  fv <- sysml_create_feature_value()
  fv_id <- fv[["@id"]]
  
  lv[["owner"]][["@id"]] <- au_id
  lv[["owningRelationship"]][["@id"]] <- fv_id
  
  fv[["memberElement"]] <- lv_id
  fv[["source"]][[1]][["@id"]] <- au_id
  fv[["target"]][[1]][["@id"]] <- lv_id
  fv[["ownedRelatedElement"]][[1]][["@id"]] <- lv_id
  fv[["owningRelatedElement"]][["@id"]] <- au_id
  
  ds[[lv_id]] <- lv
  ds[[fv_id]] <- fv
  
  ds
}

sysml_validate_mass_props_table <- function(tree, ds) {
  validate_mass_props_table(tree, ds, get_ids = sysml_get_ids, get = sysml_get_mass_props)
}

sysml_validate_mass_props_and_unc_table <- function(tree, ds) {
  validate_mass_props_and_unc_table(tree, ds, get_ids = sysml_get_ids, get = sysml_get_mass_props_and_unc)
}

sysml_set_poi_conv_from_target <- function(ds, target, mp) {
  set_poi_conv_from_target(ds, target, mp, get_by_id = sysml_get_by_id)
}

sysml_set_mass_props <- function(ds, id, mp) {
  set_mass_props(ds, id, mp, set_by_id = sysml_set_by_id)
}

sysml_set_mass_props_and_unc <- function(ds, id, mp) {
  set_mass_props_and_unc(ds, id, mp, set_by_id = sysml_set_by_id)
}

sysml_update_mass_props <- function(ds, target, sources, ...) {
  update_mass_props(ds, target, sources, set = sysml_set_mass_props, get = sysml_get_mass_props,   override = sysml_set_poi_conv_from_target
, ...)
}

sysml_update_mass_props_and_unc <- function(ds, target, sources, ...) {
  update_mass_props_and_unc(ds, target, sources, set = sysml_set_mass_props_and_unc, get = sysml_get_mass_props_and_unc, override = sysml_set_poi_conv_from_target
, ...)
}

# Test

args <- commandArgs(trailingOnly <- TRUE)
args <- c(
  "Python/trial-systemsmodeling.com/HPdK_MassPropertiesModelSmall_dump.json",
  "Python/trial-systemsmodeling.com/HPdK_MassPropertiesModelSmall_dump_rollup.json"
)

in_file <- args[1]
out_file <- args[2]

# parse json input and stack into a data frame

sysml <- read_json(in_file)
names(sysml) <- Map(f = function(i) i[["@id"]], sysml)


tree <- igraph::graph_from_edgelist(sysml_get_edgelist(sysml))

sysml_rollup <- rollupTree::rollup(
  tree,
  sysml,
  update = sysml_update_mass_props_and_unc,
  validate_ds = sysml_validate_mass_props_and_unc_table
)

write_json(unname(sysml_rollup), out_file, auto_unbox = TRUE, pretty = TRUE)

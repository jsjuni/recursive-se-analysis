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
      f = function(i) i$`@type` == "PartUsage",
      ds
    )
  ) |> unlist() |> unname()
}

sysml_get_part_with_id <- function(ds, id) {
  Filter(
    f = function(i) i$`@type` == "PartUsage" && i[["declaredShortName"]] == id,
    ds
  )[[1]]
}

sysml_get_relations_with_source <- function(ds, type, source) {
  Filter(
    f = function(i) i$`@type` == type && i[["source"]][[1]]$`@id` == source,
    ds
  )
}

sysml_get_relations_with_target <- function(ds, type, target) {
  Filter(
    f = function(i) i$`@type` == type && i[["target"]][[1]]$`@id` == target,
    ds
  )
}

sysml_get_attribute_usages_for_part <- function(ds, part, attribute) {
  Filter(
    f = function(i) i$`@type` == "AttributeUsage" && i[["declaredName"]] == attribute,
    Map(
      f = function(i) ds[[i[["target"]][[1]]$`@id`]],
      sysml_get_relations_with_source(ds, "FeatureMembership", part$`@id`)
    )
  )
}

sysml_get_negated_rational <- function(ds, operator_expression) {
  ft <- Map(
    f = function(i) ds[[unlist(i[["target"]])]],
    sysml_get_relations_with_source(ds, "ParameterMembership", operator_expression$`@id`)
  )[[1]]
  lr <- Map(
    f = function(i) ds[[unlist(i[["target"]])]],
    sysml_get_relations_with_source(ds, "FeatureValue", ft$`@id`)
  )[[1]]
  -lr[["value"]]
}

sysml_get_attribute_values <- function(ds, attribute_usage) {
  Map(
    f = function(i) {
      fv = ds[[i[["target"]][[1]]$`@id`]]
      if (fv$`@type` == "OperatorExpression") sysml_get_negated_rational(ds, fv) else fv[["value"]]
    },
    sysml_get_relations_with_source(ds, "FeatureValue", attribute_usage$`@id`)
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
        child = ds[[i[["target"]][[1]]$`@id`]][["declaredShortName"]],
        parent = ds[[i[["source"]][[1]]$`@id`]][["declaredShortName"]]
      ),
      Filter(
        f = function(i) {
          i$`@type` == "FeatureMembership" &&
            ds[[i[["source"]][[1]]$`@id`]]$`@type` == "PartUsage" &&
            ds[[i[["target"]][[1]]$`@id`]]$`@type` == "PartUsage"
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
  lr$`@type` = "LiteralRational"
  lr
}

sysml_create_literal_boolean <- function(value) {
  lr <- sysml_create_literal(value)
  lr$`@type` = "LiteralBoolean"
  lr
}

sysml_create_literal_string <- function(value) {
  lr <- sysml_create_literal(value)
  lr$`@type` = "LiteralString"
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
      `@id` = NA # Feature id
    ),
    `source` = list(
      list(
        `@id` = NA # OperatorExpression id
      )
    ),
    `target` = list(
      list(
        `@id` = NA # Feature id
      )
    ),
    `owningRelatedElement` = list(
      `@id` = NA # OperatorExpression id
    ),
    `ownedRelatedElement` = list(
      list(
        `@id` = NA # Feature id
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
    `isConstant` = FALSE,
    `isDerived` = FALSE,
    `isImpliedIncluded` = FALSE,
    `isAbstract` = FALSE,
    `isComposite` = FALSE,
    `ownedRelationship` = list(),
    `aliasIds` = list(),
    `isSufficient` = FALSE,
    `isOrdered` = FALSE,
    `elementId` = id,
    `isEnd` = FALSE,
    `isUnique` = TRUE,
    `isVariable` = FALSE,
    `owner` = list(
      `@id` = NA # OperatorExpression id
    ),
    `isPortion` = FALSE,
    `owningRelationship` = list(
      `@id` = NA # ParameterMembership or ReturnParameterMembership id
    ),
    `isLibraryElement` = FALSE
  )
}

sysml_get_owned_related_element_ids <- function(ds, id) {
  Reduce(
    f = function(l, i) union(append(l, i), sysml_get_owned_related_element_ids(ds, i)),
    x = Map(
      f = function(e) e$`@id`,
      union(
        ds[[id]][["ownedRelationship"]],
        ds[[id]][["ownedRelatedElement"]]
      )
    ),
    init = list()
  )
}

sysml_set_by_id <- function(ds, id, property, value) {
  pa <- sysml_get_part_with_id(ds, id)
  if (is.null(pa)) stop(sprintf("no part for id %s", id))
  
  au <- sysml_get_attribute_usages_for_part(ds, pa, property)[[1]]
  if (is.null(au)) stop(sprintf("no property %s for part %s", property, id))
  au_id <- au$`@id`
  
  # delete any existing values
  
  for (fv in sysml_get_relations_with_source(ds, "FeatureValue", au_id)) {
    id <- fv$`@id`
    owned <- sysml_get_owned_related_element_ids(ds, id)
    for (o in owned) {
      ds[[o]] <- NULL
    }
    ds[[id]] <- NULL
  }
  
  fv1 <- sysml_create_feature_value()
  fv1_id <- fv1$`@id`
  
  lv_or_oe  <- if (is.numeric(value)) {
    if (value < 0) {
      
      oe <- sysml_create_operator_expression()
      oe_id <- oe$`@id`
      
      pm <- sysml_create_parameter_membership()
      pm_id <- pm$`@id`

      ft1 <- sysml_create_feature()
      ft1_id <- ft1$`@id`

      fv2 <- sysml_create_feature_value()
      fv2_id <- fv2$`@id`
      
      lv <- sysml_create_literal_rational(abs(value))
      lv_id <- lv$`@id`
      
      rp <- sysml_create_return_parameter_membership()
      rp_id <- rp$`@id`
      
      ft2 <- sysml_create_feature()
      ft2_id <- ft2$`@id`
      
      oe$ownedRelationship[[1]]$`@id` <- pm_id
      oe$ownedRelationship[[1]]$`@id` <- pm_id
      oe[["ownedRelationship"]][[2]]$`@id` <- rp_id
      ds[[oe_id]] <- oe
      
      pm[["memberElement"]]$`@id` <- ft1_id
      pm[["source"]][[1]]$`@id` <- oe_id
      pm[["target"]][[1]]$`@id` <- ft1_id
      pm[["ownedRelatedElement"]][[1]]$`@id` <- ft1_id
      pm[["owningRelatedElement"]]$`@id` <- oe_id
      ds[[pm_id]] <- pm
      
      ft1[["owner"]]$`@id` <- oe_id
      ft1[["ownedRelationship"]] <- list(`@id` = fv1_id)
      ft1[["owningRelationship"]]$`@id` <- pm_id
      ds[[ft1_id]] <- ft1
      
      fv2[["memberElement"]]$`@id` <- lv_id
      fv2[["source"]][[1]]$`@id` <- ft1_id
      fv2[["target"]][[1]]$`@id` <- lv_id
      fv2[["ownedRelatedElement"]][[1]]$`@id` <- lv_id
      fv2[["owningRelatedElement"]]$`@id` <- ft1_id
      ds[[fv2_id]] <- fv2
      
      lv[["owner"]] <- ft1_id
      lv[["owningRelationship"]]$`@id` <- fv2_id
      ds[[lv_id]] <- lv
      
      rp[["memberElement"]]$`@id` <- ft2_id
      rp[["source"]][[1]]$`@id` <- oe_id
      rp[["target"]][[1]]$`@id` <- ft2_id
      rp[["owningRelatedElement"]]$`@id` <- oe_id
      rp[["ownedRelatedElement"]][[1]]$`@id` <- ft2_id
      ds[[rp_id]] <- rp

      ft2[["owner"]]$`@id` <- oe_id
      ft2[["owningRelationship"]]$`@id` <- rp_id
      ft2[["ownedRelationship"]] <- list()
      ds[[ft2_id]] <- ft2
      
      oe
       
    } else {

      lv <- sysml_create_literal_rational(value)
      lv_id <- lv$`@id`
      
      ds[[lv_id]] <- lv
       
      lv
    }
  } else if (is.logical(value)) {
    sysml_create_literal_boolean(value)
  } else if (is.character(value)) {
    sysml_create_literal_string(value)
  }
  
  lv_or_oe_id <- lv_or_oe$`@id`
  
  fv1[["memberElement"]] <- lv_or_oe_id
  fv1[["source"]][[1]]$`@id` <- au_id
  fv1[["target"]][[1]]$`@id` <- lv_or_oe_id
  fv1[["ownedRelatedElement"]][[1]]$`@id` <- lv_or_oe_id
  fv1[["owningRelatedElement"]]$`@id` <- au_id
  ds[[fv1_id]] <- fv1
  
  lv_or_oe[["owner"]]$`@id` <- au_id
  lv_or_oe[["owningRelationship"]]$`@id` <- fv1_id
  ds[[lv_or_oe_id]] <- lv_or_oe

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
names(sysml) <- Map(f = function(i) i$`@id`, sysml)


tree <- igraph::graph_from_edgelist(sysml_get_edgelist(sysml))

sysml_rollup <- rollupTree::rollup(
  tree,
  sysml,
  update = sysml_update_mass_props_and_unc,
  validate_ds = sysml_validate_mass_props_and_unc_table
)

write_json(unname(sysml_rollup), out_file, auto_unbox = TRUE, pretty = TRUE)

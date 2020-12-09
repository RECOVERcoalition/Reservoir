SELECT *
FROM molecule_hierarchy
WHERE molregno<>parent_molregno
  OR molregno <> active_molregno;

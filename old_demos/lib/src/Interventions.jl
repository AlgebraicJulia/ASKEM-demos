module Interventions
export caseload_triggered_lockdown, time_triggered_lockdown

using ModelingToolkit

function caseload_triggered_lockdown(model, cases_expression, threshold, lockdown_effect)
  ODESystem(
    [],
    ModelingToolkit.get_iv(model);
    continuous_events = [cases_expression ~ threshold] => [model.inf ~ lockdown_effect * model.inf],
    name = :caseload_triggered_lockdown
  )
end

function time_triggered_lockdown(model, T, lockdown_effect)
  t = ModelingToolkit.get_iv(model)
  ODESystem(
    [],
    t;
    discrete_events = [T] => [model.inf ~ lockdown_effect * model.inf],
    name = :time_triggered_lockdown
  )
end

end

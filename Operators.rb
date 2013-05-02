
module ConvolutionGenerator
  OPERATORS = {
    "===" => {
      :arity => 2,
      :symbol => "=",
      :languages => {
        FORTRAN => {
          :symbol => ":="
        }
      }
    }
  }
end

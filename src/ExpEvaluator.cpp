#include "ExpEvaluator.h"


/**
 * @brief Constructor initializes array pointers to NULL.
 * @details The constructor sets the interpolation scheme as the default
 *          for computing exponentials.
 */
ExpEvaluator::ExpEvaluator() {
  _interpolate = true;
  _linear_source = false;
  _exp_table = NULL;
  _polar_quad = NULL;
  _max_optical_length = MAX_OPTICAL_LENGTH;
  _exp_precision = EXP_PRECISION;
  _three_times_num_exp = 3;
}


/**
 * @brief Destructor deletes table for linear interpolation of exponentials
 */
ExpEvaluator::~ExpEvaluator() {
  if (_exp_table != NULL)
    delete [] _exp_table;
}


/**
 * @brief Set the PolarQuad to use when computing exponentials.
 * @param polar_quad a PolarQuad object pointer
 */
void ExpEvaluator::setPolarQuadrature(PolarQuad* polar_quad) {
  _polar_quad = polar_quad;
  _num_polar = _polar_quad->getNumPolarAngles();
}


/**
 * @brief Sets the maximum optical length covered in the exponential
 *        interpolation table.
 * @param max_optical_length the maximum optical length
 */
void ExpEvaluator::setMaxOpticalLength(FP_PRECISION max_optical_length) {

  if (max_optical_length <= 0)
    log_printf(ERROR, "Cannot set max optical length to %f because it "
               "must be positive.", max_optical_length);

  _max_optical_length = max_optical_length;
}


/**
 * @brief Sets the maximum acceptable approximation error for exponentials.
 * @details This routine only affects the construction of the linear
 *          interpolation table for exponentials, if in use. By default,
 *          a value of 1E-5 is used for the table, as recommended by the
 *          analysis of Yamamoto in his 2004 paper on the subject.
 * @param exp_precision the maximum exponential approximation error
 */
void ExpEvaluator::setExpPrecision(FP_PRECISION exp_precision) {

  if (exp_precision <= 0)
    log_printf(ERROR, "Cannot set exp precision to %f because it "
               "must be positive.", exp_precision);

  _exp_precision = exp_precision;
}


/**
 * @brief Use linear interpolation to compute exponentials.
 */
void ExpEvaluator::useInterpolation() {
  _interpolate = true;
}


/**
 * @brief Use the exponential intrinsic exp(...) to compute exponentials.
 */
void ExpEvaluator::useIntrinsic() {
  _interpolate = false;
}


/**
 * @brief Use linear source exponentials.
 */
void ExpEvaluator::useLinearSource() {
  _linear_source = true;
  _three_times_num_exp = 9;
}


/**
 * @brief Use flat source exponentials.
 */
void ExpEvaluator::useFlatSource() {
  _linear_source = false;
  _three_times_num_exp = 3;
}


/**
 * @brief Gets the maximum optical length covered with the exponential
 *        interpolation table.
 * @return max_optical_length the maximum optical length
 */
FP_PRECISION ExpEvaluator::getMaxOpticalLength() {
  return _max_optical_length;
}


/**
 * @brief Gets the maximum acceptable approximation error for exponentials.
 * @return the maximum exponential approximation error
 */
FP_PRECISION ExpEvaluator::getExpPrecision() {
  return _exp_precision;
}


/**
 * @brief Returns true if using linear interpolation to compute exponentials.
 * @return true if so, false otherwise
 */
bool ExpEvaluator::isUsingInterpolation() {
  return _interpolate;
}


/**
 * @brief Returns the exponential table spacing.
 * @return exponential table spacing
 */
FP_PRECISION ExpEvaluator::getTableSpacing() {
  return _exp_table_spacing;
}


FP_PRECISION ExpEvaluator::getInverseTableSpacing() {
  return _inverse_exp_table_spacing;
}


/**
 * @brief Get the number of entries in the exponential interpolation table.
 * @param entries in the interpolation table
 *
 */
int ExpEvaluator::getTableSize() {

  if (_exp_table == NULL)
    log_printf(ERROR, "Unable to return exponential table size "
               "since it has not yet been initialized");

  return _table_size;
}


/**
 * @brief Returns a pointer to the exponential interpolation table.
 * @return pointer to the exponential interpolation table
 */
FP_PRECISION* ExpEvaluator::getExpTable() {

  if (_exp_table == NULL)
    log_printf(ERROR, "Unable to return exponential table "
               "since it has not yet been initialized");

  return _exp_table;
}


/**
 * @brief If using linear interpolation, builds the table for each polar angle.
 * @param tolerance the minimum acceptable interpolation accuracy
 */
void ExpEvaluator::initialize() {

  log_printf(INFO, "Initializing exponential evaluator...");

  _sin_theta = _polar_quad->getSinThetas();
  _inv_sin_theta = _polar_quad->getInverseSinThetas();

  /* Set size of interpolation table */
  int num_array_values;
  if (_linear_source)
    num_array_values = _max_optical_length * sqrt(1. / (8. * _exp_precision));
  else
    num_array_values = _max_optical_length *
        pow(1. / (72. * sqrt(3.0) * _exp_precision), 1.0/3.0);

  _exp_table_spacing = _max_optical_length / num_array_values;

  /* Increment the number of vaues in the array to ensure that a tau equal to
   * max_optical_length resides as the final entry in the table */
  num_array_values += 1;

  /* Compute the reciprocal of the table entry spacing */
  _inverse_exp_table_spacing = 1.0 / _exp_table_spacing;

  /* Allocate array for the table */
  if (_exp_table != NULL)
    delete [] _exp_table;

  _table_size = _three_times_num_exp * _num_polar * num_array_values;
  _exp_table = new FP_PRECISION[_table_size];

  FP_PRECISION expon;
  FP_PRECISION intercept;
  FP_PRECISION slope_F1;
  FP_PRECISION st, ist, ist2, ist3;
  FP_PRECISION tau_a, tau_m, tau_a2, tau_a3;
  FP_PRECISION itau_a, itau_a2, itau_a3;
  FP_PRECISION exp_const_1, exp_const_2, exp_const_3;
  int index;

  /* Create exponential linear interpolation table */
  for (int i=0; i < num_array_values; i++) {
    for (int p=0; p < _num_polar; p++) {

      index = _three_times_num_exp * (_num_polar * i + p);

      st = _sin_theta[p];
      ist = _inv_sin_theta[p];
      ist2 = ist * ist;
      ist3 = ist * ist * ist;

      tau_a = i * _exp_table_spacing;
      tau_a2 = tau_a * tau_a;
      tau_a3 = tau_a * tau_a * tau_a;
      itau_a = 1.0 / tau_a;
      itau_a2 = itau_a * itau_a;
      itau_a3 = itau_a * itau_a * itau_a;
      tau_m = tau_a * ist;
      expon = exp(- tau_m);

      /* Compute F1 */
      exp_const_1 = (1.0 - expon);
      exp_const_2 = expon * ist;
      exp_const_3 = -0.5 * exp_const_2 * ist;
      _exp_table[index    ] = exp_const_1;
      _exp_table[index + 1] = exp_const_2;
      _exp_table[index + 2] = exp_const_3;

      if (_linear_source) {

        /* Compute F2 */
        if (tau_a == 0.0) {
          exp_const_1 = 0.0;
          exp_const_2 = 0.0;
          exp_const_3 = ist3 / 6.0;
        }
        else {
          exp_const_1 = 2.0 * itau_a * (expon - 1) + ist * (expon + 1);
          exp_const_2 = 2.0 * itau_a2 +
              expon * (- ist2 - 2.0 * ist * itau_a - 2.0 * itau_a2);
          exp_const_3 = 2.0 * itau_a3 +
              expon * (- ist3 / 2.0 - ist2 * itau_a - 2.0 * itau_a2 * ist -
                       2.0 * itau_a3);
        }
        _exp_table[index + 3] = exp_const_1;
        _exp_table[index + 4] = exp_const_2;
        _exp_table[index + 5] = exp_const_3;

        /* Compute H */
        if (tau_a == 0.0) {
          exp_const_1 = 0.0;
          exp_const_2 = 0.5 * ist;
          exp_const_3 = -1.0 * ist * ist / 3.0;
        }
        else {
          exp_const_1 = (-expon * (tau_a + st) + st) / tau_a;
          exp_const_2 = (expon * (tau_a * tau_a + tau_a * st + st * st)
                         - st * st) / (tau_a * tau_a * st);
          exp_const_3 = 1.0 / (2 * tau_a * tau_a * tau_a * st * st) *
              (-expon * (tau_a * tau_a * tau_a + tau_a * tau_a * st +
                         2 * tau_a * st * st + 2 * st * st * st) +
               2 * st * st * st);
        }
        _exp_table[index + 6] = exp_const_1;
        _exp_table[index + 7] = exp_const_2;
        _exp_table[index + 8] = exp_const_3;
      }
    }
  }
}


/**
 * @brief Computes the G2 exponential term for a optical length and polar angle.
 * @details This method computes the H exponential term from Ferrer [1]
 *          for some optical path length and polar angle. This method
 *          uses either a linear interpolation table (default) or the
 *          exponential intrinsic exp(...) function.
 *
 *            [1] R. Ferrer and J. Rhodes III, "A Linear Source Approximation
 *                Scheme for the Method of Characteristics", Nuclear Science and
 *                Engineering, Volume 182, February 2016.
 *
 * @param tau the optical path length (e.g., sigma_t times length)
 * @param polar the polar angle index
 * @return the evaluated exponential
 */
FP_PRECISION ExpEvaluator::computeExponentialG2(FP_PRECISION tau, int polar) {

  tau = std::max(tau, FP_PRECISION(1.e-10));

  if (tau == 0.0)
    return 0.0;

  FP_PRECISION tau_m = tau * _inv_sin_theta[polar];
  return _inv_sin_theta[polar] / 3.0 - (0.5 / tau + 1.0 / (tau_m * tau))
      * (1.0 + tau_m / 2.0 - (1.0 + 1.0 / tau_m) *
         (1.0 - exp(- tau_m)));
}

//// Type: UnityEngine.Mathd
//// Assembly: UnityEngine, Version=0.0.0.0, Culture=neutral, PublicKeyToken=null
//// Assembly location: C:\Program Files (x86)\Unity\Editor\Data\Managed\UnityEngine.dll
//using System;
//using System.Runtime.CompilerServices;
//
//
//namespace UnityEngine {
//    public struct MathDec {
//        public const decimal PI = 3.141593m;
//        public const decimal Deg2Rad = 0.01745329m;
//        public const decimal Rad2Deg = 57.29578m;
//        public const decimal Epsilon = 1.401298E-45m;
//
//        public static decimal Sin(decimal d) {
//            return Math.Sin(d);
//        }
//
//        public static decimal Cos(decimal d) {
//            return Math.Cos(d);
//        }
//
//        public static decimal Tan(decimal d) {
//            return Math.Tan(d);
//        }
//
//        public static decimal Asin(decimal d) {
//            return Math.Asin(d);
//        }
//
//        public static decimal Acos(decimal d) {
//            return Math.Acos(d);
//        }
//
//        public static decimal Atan(decimal d) {
//            return Math.Atan(d);
//        }
//
//        public static decimal Atan2(decimal y, decimal x) {
//            return Math.Atan2(y, x);
//        }
//
//        public static decimal Sqrt(decimal d) {
//            return Math.Sqrt(d);
//        }
//
//        public static decimal Abs(decimal d) {
//            return Math.Abs(d);
//        }
//
//        public static int Abs(int value) {
//            return Math.Abs(value);
//        }
//
//        public static decimal Min(decimal a, decimal b) {
//            if (a < b)
//                return a;
//            else
//                return b;
//        }
//
//        public static decimal Min(params decimal[] values) {
//            int length = values.Length;
//            if (length == 0)
//                return 0.0d;
//            decimal num = values[0];
//            for (int index = 1; index < length; ++index) {
//                if (values[index] < num)
//                    num = values[index];
//            }
//            return num;
//        }
//
//        public static int Min(int a, int b) {
//            if (a < b)
//                return a;
//            else
//                return b;
//        }
//
//        public static int Min(params int[] values) {
//            int length = values.Length;
//            if (length == 0)
//                return 0;
//            int num = values[0];
//            for (int index = 1; index < length; ++index) {
//                if (values[index] < num)
//                    num = values[index];
//            }
//            return num;
//        }
//
//        public static decimal Max(decimal a, decimal b) {
//            if (a > b)
//                return a;
//            else
//                return b;
//        }
//
//        public static decimal Max(params decimal[] values) {
//            int length = values.Length;
//            if (length == 0)
//                return 0d;
//            decimal num = values[0];
//            for (int index = 1; index < length; ++index) {
//                if ((decimal)values[index] > (decimal)num)
//                    num = values[index];
//            }
//            return num;
//        }
//
//        public static int Max(int a, int b) {
//            if (a > b)
//                return a;
//            else
//                return b;
//        }
//
//        public static int Max(params int[] values) {
//            int length = values.Length;
//            if (length == 0)
//                return 0;
//            int num = values[0];
//            for (int index = 1; index < length; ++index) {
//                if (values[index] > num)
//                    num = values[index];
//            }
//            return num;
//        }
//
//        public static decimal Pow(decimal d, decimal p) {
//            return Math.Pow(d, p);
//        }
//
//        public static decimal Exp(decimal power) {
//            return Math.Exp(power);
//        }
//
//        public static decimal Log(decimal d, decimal p) {
//            return Math.Log(d, p);
//        }
//
//        public static decimal Log(decimal d) {
//            return Math.Log(d);
//        }
//
//        public static decimal Log10(decimal d) {
//            return Math.Log10(d);
//        }
//
//        public static decimal Ceil(decimal d) {
//            return Math.Ceiling(d);
//        }
//
//        public static decimal Floor(decimal d) {
//            return Math.Floor(d);
//        }
//
//        public static decimal Round(decimal d) {
//            return Math.Round(d);
//        }
//
//        public static int CeilToInt(decimal d) {
//            return (int)Math.Ceiling(d);
//        }
//
//        public static int FloorToInt(decimal d) {
//            return (int)Math.Floor(d);
//        }
//
//        public static int RoundToInt(decimal d) {
//            return (int)Math.Round(d);
//        }
//
//        public static decimal Sign(decimal d) {
//            return d >= 0.0 ? 1d : -1d;
//        }
//
//        public static decimal Clamp(decimal value, decimal min, decimal max) {
//            if (value < min)
//                value = min;
//            else if (value > max)
//                value = max;
//            return value;
//        }
//
//        public static int Clamp(int value, int min, int max) {
//            if (value < min)
//                value = min;
//            else if (value > max)
//                value = max;
//            return value;
//        }
//
//        public static decimal Clamp01(decimal value) {
//            if (value < 0.0)
//                return 0.0d;
//            if (value > 1.0)
//                return 1d;
//            else
//                return value;
//        }
//
//        public static decimal Lerp(decimal from, decimal to, decimal t) {
//            return from + (to - from) * Mathd.Clamp01(t);
//        }
//
//        public static decimal LerpAngle(decimal a, decimal b, decimal t) {
//            decimal num = Mathd.Repeat(b - a, 360d);
//            if (num > 180.0d)
//                num -= 360d;
//            return a + num * Mathd.Clamp01(t);
//        }
//
//        public static decimal MoveTowards(decimal current, decimal target, decimal maxDelta) {
//            if (Mathd.Abs(target - current) <= maxDelta)
//                return target;
//            else
//                return current + Mathd.Sign(target - current) * maxDelta;
//        }
//
//        public static decimal MoveTowardsAngle(decimal current, decimal target, decimal maxDelta) {
//            target = current + Mathd.DeltaAngle(current, target);
//            return Mathd.MoveTowards(current, target, maxDelta);
//        }
//
//        public static decimal SmoothStep(decimal from, decimal to, decimal t) {
//            t = Mathd.Clamp01(t);
//            t = (-2.0 * t * t * t + 3.0 * t * t);
//            return to * t + from * (1.0 - t);
//        }
//
//        public static decimal Gamma(decimal value, decimal absmax, decimal gamma) {
//            bool flag = false;
//            if (value < 0.0)
//                flag = true;
//            decimal num1 = Mathd.Abs(value);
//            if (num1 > absmax) {
//                if (flag)
//                    return -num1;
//                else
//                    return num1;
//            } else {
//                decimal num2 = Mathd.Pow(num1 / absmax, gamma) * absmax;
//                if (flag)
//                    return -num2;
//                else
//                    return num2;
//            }
//        }
//
//        public static bool Approximately(decimal a, decimal b) {
//            return Mathd.Abs(b - a) < Mathd.Max(1E-06d * Mathd.Max(Mathd.Abs(a), Mathd.Abs(b)), 1.121039E-44d);
//        }
//		/**
//		 * In this version, the tolerance is based on a fixed precision. 
//		 * For example, the values may be specified in units of degrees, or millimeters, etc.
//		 * Use this version when you know how precise you want to be.
//		 */
//		public static bool ApproximatelyFixed(decimal a, decimal b, decimal baseValue) {
//			decimal tolerance = baseValue * DoublePrecision;
//			decimal difference = Mathd.Abs(a - b);
//			return difference <= tolerance;
//		}
//		/**
//		 * In this version, the tolerance is based on the values themselves.
//		 * This version is used when the caller is blind with regards to the meaning of the values
//		 * and the degree of necessary precision. The toleranceMultiplier is more or less how many
//		 * "steps" the values are allowed to be apart. toleranceMultipliers less than 1 translate essentially
//		 * to "must be exactly equal", and values less than 1 will always fail.
//		 */
//		public static bool ApproximatelyRelative(decimal a, decimal b, decimal toleranceMultiplier) {
//			decimal tolerance = toleranceMultiplier * DoublePrecision * Mathd.Max(Mathd.Abs(a), Mathd.Abs(b));
//			decimal difference = Mathd.Abs(a - b);
//			return difference <= tolerance;
//		}
//
//        public static decimal SmoothDamp(decimal current, decimal target, ref decimal currentVelocity, decimal smoothTime, decimal maxSpeed) {
//            decimal deltaTime = (decimal)Time.deltaTime;
//            return Mathd.SmoothDamp(current, target, ref currentVelocity, smoothTime, maxSpeed, deltaTime);
//        }
//
//        public static decimal SmoothDamp(decimal current, decimal target, ref decimal currentVelocity, decimal smoothTime) {
//            decimal deltaTime = Time.deltaTime;
//            decimal maxSpeed = decimal.PositiveInfinity;
//            return Mathd.SmoothDamp(current, target, ref currentVelocity, smoothTime, maxSpeed, deltaTime);
//        }
//
//        public static decimal SmoothDamp(decimal current, decimal target, ref decimal currentVelocity, decimal smoothTime, decimal maxSpeed, decimal deltaTime) {
//            smoothTime = Mathd.Max(0.0001d, smoothTime);
//            decimal num1 = 2d / smoothTime;
//            decimal num2 = num1 * deltaTime;
//            decimal num3 = (1.0d / (1.0d + num2 + 0.479999989271164d * num2 * num2 + 0.234999999403954d * num2 * num2 * num2));
//            decimal num4 = current - target;
//            decimal num5 = target;
//            decimal max = maxSpeed * smoothTime;
//            decimal num6 = Mathd.Clamp(num4, -max, max);
//            target = current - num6;
//            decimal num7 = (currentVelocity + num1 * num6) * deltaTime;
//            currentVelocity = (currentVelocity - num1 * num7) * num3;
//            decimal num8 = target + (num6 + num7) * num3;
//            if (num5 - current > 0.0 == num8 > num5) {
//                num8 = num5;
//                currentVelocity = (num8 - num5) / deltaTime;
//            }
//            return num8;
//        }
//
//        public static decimal SmoothDampAngle(decimal current, decimal target, ref decimal currentVelocity, decimal smoothTime, decimal maxSpeed) {
//            decimal deltaTime = (decimal)Time.deltaTime;
//            return Mathd.SmoothDampAngle(current, target, ref currentVelocity, smoothTime, maxSpeed, deltaTime);
//        }
//
//        public static decimal SmoothDampAngle(decimal current, decimal target, ref decimal currentVelocity, decimal smoothTime) {
//            decimal deltaTime = (decimal)Time.deltaTime;
//            decimal maxSpeed = decimal.PositiveInfinity;
//            return Mathd.SmoothDampAngle(current, target, ref currentVelocity, smoothTime, maxSpeed, deltaTime);
//        }
//
//        public static decimal SmoothDampAngle(decimal current, decimal target, ref decimal currentVelocity, decimal smoothTime, decimal maxSpeed, decimal deltaTime) {
//            target = current + Mathd.DeltaAngle(current, target);
//            return Mathd.SmoothDamp(current, target, ref currentVelocity, smoothTime, maxSpeed, deltaTime);
//        }
//
//        public static decimal Repeat(decimal t, decimal length) {
//            return t - Mathd.Floor(t / length) * length;
//        }
//
//        public static decimal PingPong(decimal t, decimal length) {
//            t = Mathd.Repeat(t, length * 2d);
//            return length - Mathd.Abs(t - length);
//        }
//
//        public static decimal InverseLerp(decimal from, decimal to, decimal value) {
//            if (from < to) {
//                if (value < from)
//                    return 0d;
//                if (value > to)
//                    return 1d;
//                value -= from;
//                value /= to - from;
//                return value;
//            } else {
//                if (from <= to)
//                    return 0d;
//                if (value < to)
//                    return 1d;
//                if (value > from)
//                    return 0d;
//                else
//                    return (1.0d - (value - to) / (from - to));
//            }
//        }
//
//        public static decimal DeltaAngle(decimal current, decimal target) {
//            decimal num = Mathd.Repeat(target - current, 360d);
//            if (num > 180.0d)
//                num -= 360d;
//            return num;
//        }
//
//        internal static bool LineIntersection(Vector2d p1, Vector2d p2, Vector2d p3, Vector2d p4, ref Vector2d result) {
//            decimal num1 = p2.x - p1.x;
//            decimal num2 = p2.y - p1.y;
//            decimal num3 = p4.x - p3.x;
//            decimal num4 = p4.y - p3.y;
//            decimal num5 = num1 * num4 - num2 * num3;
//            if (num5 == 0.0d)
//                return false;
//            decimal num6 = p3.x - p1.x;
//            decimal num7 = p3.y - p1.y;
//            decimal num8 = (num6 * num4 - num7 * num3) / num5;
//            result = new Vector2d(p1.x + num8 * num1, p1.y + num8 * num2);
//            return true;
//        }
//
//        internal static bool LineSegmentIntersection(Vector2d p1, Vector2d p2, Vector2d p3, Vector2d p4, ref Vector2d result) {
//            decimal num1 = p2.x - p1.x;
//            decimal num2 = p2.y - p1.y;
//            decimal num3 = p4.x - p3.x;
//            decimal num4 = p4.y - p3.y;
//            decimal num5 = (num1 * num4 - num2 * num3);
//            if (num5 == 0.0d)
//                return false;
//            decimal num6 = p3.x - p1.x;
//            decimal num7 = p3.y - p1.y;
//            decimal num8 = (num6 * num4 - num7 * num3) / num5;
//            if (num8 < 0.0d || num8 > 1.0d)
//                return false;
//            decimal num9 = (num6 * num2 - num7 * num1) / num5;
//            if (num9 < 0.0d || num9 > 1.0d)
//                return false;
//            result = new Vector2d(p1.x + num8 * num1, p1.y + num8 * num2);
//            return true;
//        }
//    }
//}
//

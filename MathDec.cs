// Type: UnityEngine.MathDec
// Assembly: UnityEngine, Version=0.0.0.0, Culture=neutral, PublicKeyToken=null
// Assembly location: C:\Program Files (x86)\Unity\Editor\Data\Managed\UnityEngine.dll
using System;
using System.Runtime.CompilerServices;


namespace UnityEngine {
    public struct MathDec {
		public const decimal PI = 3.1415926535897932384626433832795028841972m;
		public const decimal Deg2Rad = 0.0174532925199432957692369076848861271344m;
		public const decimal Rad2Deg = 57.2957795130823208767981548141051703324049m;
        public const decimal Epsilon = 1.401298E-45m;
		public const decimal DefaultApproximatelyValue = 1.121039E-44m;
		public const decimal DoublePrecision = 2.22044604925031E-16m;	//	Minimum relative value. Anything smaller than this gets rounded off. In most cases, it might be useful to double this.
		public const decimal DoublePrecisionMultiplier = 4503599627370500m;
		public const decimal DecimalPrecision = 2.22044604925031E-16m;	//	Minimum relative value. Anything smaller than this gets rounded off. In most cases, it might be useful to double this.
		public const decimal DecimalPrecisionMultiplier = 4503599627370500m;

        public static decimal Sin(decimal d) {
			return (decimal)Math.Sin((double)d);
        }

        public static decimal Cos(decimal d) {
			return (decimal)Math.Cos((double)d);
        }

        public static decimal Tan(decimal d) {
			return (decimal)Math.Tan((double)d);
        }

        public static decimal Asin(decimal d) {
			return (decimal)Math.Asin((double)d);
        }

        public static decimal Acos(decimal d) {
			return (decimal)Math.Acos((double)d);
        }

        public static decimal Atan(decimal d) {
			return (decimal)Math.Atan((double)d);
        }

        public static decimal Atan2(decimal y, decimal x) {
			return (decimal)Math.Atan2((double)y, (double)x);
        }
		public static decimal Sqrt(decimal d, decimal epsilon){
			if (d < 0m) throw new OverflowException("Cannot calculate square root from a negative number");
			
			decimal current = (decimal)Math.Sqrt((double)d);
			decimal previous;
			do
			{
				previous = current;
				if (previous == 0.0M) return 0;
				current = (previous + d / previous) / 2;
			}
			while (Math.Abs(previous - current) > epsilon);
			return current;
		}
        public static decimal Sqrt(decimal d) {
			return Sqrt (d, Epsilon);
//			return (decimal)Math.Sqrt((double)d);
        }

        public static decimal Abs(decimal d) {
			return (decimal)Math.Abs((double)d);
        }

        public static int Abs(int value) {
            return Math.Abs(value);
        }

        public static decimal Min(decimal a, decimal b) {
            if (a < b)
                return a;
            else
                return b;
        }

        public static decimal Min(params decimal[] values) {
            int length = values.Length;
            if (length == 0)
                return 0.0m;
            decimal num = values[0];
            for (int index = 1; index < length; ++index) {
                if (values[index] < num)
                    num = values[index];
            }
            return num;
        }

        public static int Min(int a, int b) {
            if (a < b)
                return a;
            else
                return b;
        }

        public static int Min(params int[] values) {
            int length = values.Length;
            if (length == 0)
                return 0;
            int num = values[0];
            for (int index = 1; index < length; ++index) {
                if (values[index] < num)
                    num = values[index];
            }
            return num;
        }

        public static decimal Max(decimal a, decimal b) {
            if (a > b)
                return a;
            else
                return b;
        }

        public static decimal Max(params decimal[] values) {
            int length = values.Length;
            if (length == 0)
                return 0m;
            decimal num = values[0];
            for (int index = 1; index < length; ++index) {
                if ((decimal)values[index] > (decimal)num)
                    num = values[index];
            }
            return num;
        }

        public static int Max(int a, int b) {
            if (a > b)
                return a;
            else
                return b;
        }

        public static int Max(params int[] values) {
            int length = values.Length;
            if (length == 0)
                return 0;
            int num = values[0];
            for (int index = 1; index < length; ++index) {
                if (values[index] > num)
                    num = values[index];
            }
            return num;
        }

        public static decimal Pow(decimal d, decimal p) {
			return (decimal)Math.Pow((double)d, (double)p);
        }

        public static decimal Exp(decimal power) {
			return (decimal)Math.Exp((double)power);
        }

        public static decimal Log(decimal d, decimal p) {
			return (decimal)Math.Log((double)d, (double)p);
        }

        public static decimal Log(decimal d) {
			return (decimal)Math.Log((double)d);
        }

        public static decimal Log10(decimal d) {
			return (decimal)Math.Log10((double)d);
        }

        public static decimal Ceil(decimal d) {
			return Math.Ceiling(d);
        }

        public static decimal Floor(decimal d) {
			return Math.Floor(d);
        }

        public static decimal Round(decimal d) {
			return Math.Round(d);
        }

        public static int CeilToInt(decimal d) {
            return (int)Math.Ceiling(d);
        }

        public static int FloorToInt(decimal d) {
            return (int)Math.Floor(d);
        }

        public static int RoundToInt(decimal d) {
            return (int)Math.Round(d);
        }

        public static decimal Sign(decimal d) {
            return d >= 0.0m ? 1m : -1m;
        }

        public static decimal Clamp(decimal value, decimal min, decimal max) {
            if (value < min)
                value = min;
            else if (value > max)
                value = max;
            return value;
        }

        public static int Clamp(int value, int min, int max) {
            if (value < min)
                value = min;
            else if (value > max)
                value = max;
            return value;
        }

        public static decimal Clamp01(decimal value) {
            if (value < 0.0m)
                return 0.0m;
            if (value > 1.0m)
                return 1m;
            else
                return value;
        }

        public static decimal Lerp(decimal from, decimal to, decimal t) {
            return from + (to - from) * MathDec.Clamp01(t);
        }

        public static decimal LerpAngle(decimal a, decimal b, decimal t) {
            decimal num = MathDec.Repeat(b - a, 360m);
            if (num > 180.0m)
                num -= 360m;
            return a + num * MathDec.Clamp01(t);
        }

        public static decimal MoveTowards(decimal current, decimal target, decimal maxDelta) {
            if (MathDec.Abs(target - current) <= maxDelta)
                return target;
            else
                return current + MathDec.Sign(target - current) * maxDelta;
        }

        public static decimal MoveTowardsAngle(decimal current, decimal target, decimal maxDelta) {
            target = current + MathDec.DeltaAngle(current, target);
            return MathDec.MoveTowards(current, target, maxDelta);
        }

        public static decimal SmoothStep(decimal from, decimal to, decimal t) {
            t = MathDec.Clamp01(t);
            t = (-2.0m * t * t * t + 3.0m * t * t);
            return to * t + from * (1.0m - t);
        }

        public static decimal Gamma(decimal value, decimal absmax, decimal gamma) {
            bool flag = false;
            if (value < 0.0m)
                flag = true;
            decimal num1 = MathDec.Abs(value);
            if (num1 > absmax) {
                if (flag)
                    return -num1;
                else
                    return num1;
            } else {
                decimal num2 = MathDec.Pow(num1 / absmax, gamma) * absmax;
                if (flag)
                    return -num2;
                else
                    return num2;
            }
        }

        public static bool Approximately(decimal a, decimal b) {
            return MathDec.Abs(b - a) < MathDec.Max(1E-06m * MathDec.Max(MathDec.Abs(a), MathDec.Abs(b)), DefaultApproximatelyValue);
        }
		/**
		 * In this version, the tolerance is based on a fixed precision. 
		 * For example, the values may be specified in units of degrees, or millimeters, etc.
		 * Use this version when you know how precise you want to be.
		 */
		public static bool ApproximatelyFixed(decimal a, decimal b, decimal baseValue) {
			decimal tolerance = baseValue * DoublePrecision;
			decimal difference = MathDec.Abs(a - b);
			return difference <= tolerance;
		}
		/**
		 * In this version, the tolerance is based on the values themselves.
		 * This version is used when the caller is blind with regards to the meaning of the values
		 * and the degree of necessary precision. The toleranceMultiplier is more or less how many
		 * "steps" the values are allowed to be apart. toleranceMultipliers less than 1 translate essentially
		 * to "must be exactly equal", and values less than 1 will always fail.
		 */
		public static bool ApproximatelyRelative(decimal a, decimal b, decimal toleranceMultiplier) {
			decimal tolerance = toleranceMultiplier * DoublePrecision * MathDec.Max(MathDec.Abs(a), MathDec.Abs(b));
			decimal difference = MathDec.Abs(a - b);
			return difference <= tolerance;
		}

        public static decimal SmoothDamp(decimal current, decimal target, ref decimal currentVelocity, decimal smoothTime, decimal maxSpeed) {
            decimal deltaTime = (decimal)Time.deltaTime;
            return MathDec.SmoothDamp(current, target, ref currentVelocity, smoothTime, maxSpeed, deltaTime);
        }

        public static decimal SmoothDamp(decimal current, decimal target, ref decimal currentVelocity, decimal smoothTime) {
            decimal deltaTime = (decimal)Time.deltaTime;
            decimal maxSpeed = decimal.MaxValue;
            return MathDec.SmoothDamp(current, target, ref currentVelocity, smoothTime, maxSpeed, deltaTime);
        }

        public static decimal SmoothDamp(decimal current, decimal target, ref decimal currentVelocity, decimal smoothTime, decimal maxSpeed, decimal deltaTime) {
            smoothTime = MathDec.Max(0.0001m, smoothTime);
            decimal num1 = 2m / smoothTime;
            decimal num2 = num1 * deltaTime;
            decimal num3 = (1.0m / (1.0m + num2 + 0.479999989271164m * num2 * num2 + 0.234999999403954m * num2 * num2 * num2));
            decimal num4 = current - target;
            decimal num5 = target;
            decimal max = maxSpeed * smoothTime;
            decimal num6 = MathDec.Clamp(num4, -max, max);
            target = current - num6;
            decimal num7 = (currentVelocity + num1 * num6) * deltaTime;
            currentVelocity = (currentVelocity - num1 * num7) * num3;
            decimal num8 = target + (num6 + num7) * num3;
            if (num5 - current > 0.0m == num8 > num5) {
                num8 = num5;
                currentVelocity = (num8 - num5) / deltaTime;
            }
            return num8;
        }

        public static decimal SmoothDampAngle(decimal current, decimal target, ref decimal currentVelocity, decimal smoothTime, decimal maxSpeed) {
            decimal deltaTime = (decimal)Time.deltaTime;
            return MathDec.SmoothDampAngle(current, target, ref currentVelocity, smoothTime, maxSpeed, deltaTime);
        }

        public static decimal SmoothDampAngle(decimal current, decimal target, ref decimal currentVelocity, decimal smoothTime) {
            decimal deltaTime = (decimal)Time.deltaTime;
            decimal maxSpeed = decimal.MaxValue;
            return MathDec.SmoothDampAngle(current, target, ref currentVelocity, smoothTime, maxSpeed, deltaTime);
        }

        public static decimal SmoothDampAngle(decimal current, decimal target, ref decimal currentVelocity, decimal smoothTime, decimal maxSpeed, decimal deltaTime) {
            target = current + MathDec.DeltaAngle(current, target);
            return MathDec.SmoothDamp(current, target, ref currentVelocity, smoothTime, maxSpeed, deltaTime);
        }

        public static decimal Repeat(decimal t, decimal length) {
            return t - MathDec.Floor(t / length) * length;
        }

        public static decimal PingPong(decimal t, decimal length) {
            t = MathDec.Repeat(t, length * 2m);
            return length - MathDec.Abs(t - length);
        }

        public static decimal InverseLerp(decimal from, decimal to, decimal value) {
            if (from < to) {
                if (value < from)
                    return 0m;
                if (value > to)
                    return 1m;
                value -= from;
                value /= to - from;
                return value;
            } else {
                if (from <= to)
                    return 0m;
                if (value < to)
                    return 1m;
                if (value > from)
                    return 0m;
                else
                    return (1.0m - (value - to) / (from - to));
            }
        }

        public static decimal DeltaAngle(decimal current, decimal target) {
            decimal num = MathDec.Repeat(target - current, 360m);
            if (num > 180.0m)
                num -= 360m;
            return num;
        }

        internal static bool LineIntersection(DecVector2 p1, DecVector2 p2, DecVector2 p3, DecVector2 p4, ref DecVector2 result) {
            decimal num1 = p2.x - p1.x;
            decimal num2 = p2.y - p1.y;
            decimal num3 = p4.x - p3.x;
            decimal num4 = p4.y - p3.y;
            decimal num5 = num1 * num4 - num2 * num3;
            if (num5 == 0.0m)
                return false;
            decimal num6 = p3.x - p1.x;
            decimal num7 = p3.y - p1.y;
            decimal num8 = (num6 * num4 - num7 * num3) / num5;
            result = new DecVector2(p1.x + num8 * num1, p1.y + num8 * num2);
            return true;
        }

        internal static bool LineSegmentIntersection(DecVector2 p1, DecVector2 p2, DecVector2 p3, DecVector2 p4, ref DecVector2 result) {
            decimal num1 = p2.x - p1.x;
            decimal num2 = p2.y - p1.y;
            decimal num3 = p4.x - p3.x;
            decimal num4 = p4.y - p3.y;
            decimal num5 = (num1 * num4 - num2 * num3);
            if (num5 == 0.0m)
                return false;
            decimal num6 = p3.x - p1.x;
            decimal num7 = p3.y - p1.y;
            decimal num8 = (num6 * num4 - num7 * num3) / num5;
            if (num8 < 0.0m || num8 > 1.0m)
                return false;
            decimal num9 = (num6 * num2 - num7 * num1) / num5;
            if (num9 < 0.0m || num9 > 1.0m)
                return false;
            result = new DecVector2(p1.x + num8 * num1, p1.y + num8 * num2);
            return true;
        }
    }
}


// Type: UnityEngine.DecVector3
// Assembly: UnityEngine, Version=0.0.0.0, Culture=neutral, PublicKeyToken=null
// Assembly location: C:\Program Files (x86)\Unity\Editor\Data\Managed\UnityEngine.dll
using System;
using System.Runtime.CompilerServices;

namespace UnityEngine {
    public struct DecVector3 {
		public const float kEpsilon = 1E-05f;
		public const decimal degreesPerRadian = 57.2957795130823m;
		public const string defaultNumberFormat = "G5";
        public decimal x;
        public decimal y;
        public decimal z;

		public static readonly DecVector3 up = new DecVector3(0m,1m,0m);

        public decimal this[int index] {
            get {
                switch (index) {
                    case 0:
                        return this.x;
                    case 1:
                        return this.y;
                    case 2:
                        return this.z;
                    default:
                        throw new IndexOutOfRangeException("Invalid index!");
                }
            }
            set {
                switch (index) {
                    case 0:
                        this.x = value;
                        break;
                    case 1:
                        this.y = value;
                        break;
                    case 2:
                        this.z = value;
                        break;
                    default:
                        throw new IndexOutOfRangeException("Invalid DecVector3 index!");
                }
            }
        }

        public DecVector3 normalized {
            get {
                return DecVector3.Normalize(this);
            }
        }

        public decimal magnitude {
            get {
                return Math.Sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
            }
        }

        public decimal sqrMagnitude {
            get {
                return this.x * this.x + this.y * this.y + this.z * this.z;
            }
        }

        public static DecVector3 zero {
            get {
                return new DecVector3(0m, 0m, 0m);
            }
        }

        public static DecVector3 one {
            get {
                return new DecVector3(1m, 1m, 1m);
            }
        }

        public static DecVector3 forward {
            get {
                return new DecVector3(0m, 0m, 1m);
            }
        }

        public static DecVector3 back {
            get {
                return new DecVector3(0m, 0m, -1m);
            }
        }

//        public static DecVector3 up {
//            get {
//                return new DecVector3(0m, 1m, 0m);
//            }
//        }

        public static DecVector3 down {
            get {
                return new DecVector3(0m, -1m, 0m);
            }
        }

        public static DecVector3 left {
            get {
                return new DecVector3(-1m, 0m, 0m);
            }
        }

        public static DecVector3 right {
            get {
                return new DecVector3(1m, 0m, 0m);
            }
        }

        public DecVector3(decimal x, decimal y, decimal z) {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        public DecVector3(float x, float y, float z) {
            this.x = (decimal)x;
            this.y = (decimal)y;
            this.z = (decimal)z;
        }

        public DecVector3(Vector3 v3) {
            this.x = (decimal)v3.x;
            this.y = (decimal)v3.y;
            this.z = (decimal)v3.z;
        }

        public DecVector3(decimal x, decimal y) {
            this.x = x;
            this.y = y;
            this.z = 0m;
        }

        public static DecVector3 operator +(DecVector3 a, DecVector3 b) {
            return new DecVector3(a.x + b.x, a.y + b.y, a.z + b.z);
        }

        public static DecVector3 operator -(DecVector3 a, DecVector3 b) {
            return new DecVector3(a.x - b.x, a.y - b.y, a.z - b.z);
        }

        public static DecVector3 operator -(DecVector3 a) {
            return new DecVector3(-a.x, -a.y, -a.z);
        }

        public static DecVector3 operator *(DecVector3 a, decimal d) {
            return new DecVector3(a.x * d, a.y * d, a.z * d);
        }

        public static DecVector3 operator *(decimal d, DecVector3 a) {
            return new DecVector3(a.x * d, a.y * d, a.z * d);
        }

        public static DecVector3 operator /(DecVector3 a, decimal d) {
            return new DecVector3(a.x / d, a.y / d, a.z / d);
        }

        public static bool operator ==(DecVector3 lhs, DecVector3 rhs) {
			//	What's up with the division operation?
			//	Wait ... less than zero? ... what? I have zero comprehension of this. Shouldn't this always fail? What exactly is this testing?
			//	Ultimately, isn't this just testing that all components are equal? (Not quite, this would eliminate really tiny fractions)
			//	Okay, one thing this would do is eliminate really tiny fractions. Except ... no, it wouldn't. I think this is just stupid.

			//	Yes ... this always fails. Even when testing an instance against itself. That should be "duh" test case.
			//	I'm betting he wrote it this way so he could adjust the precision arbitrarily. Dumb.

			return lhs.x == rhs.x 
				&& lhs.y == rhs.y 
				&& lhs.z == rhs.z;
//			Debug.Log ("'==' of " + lhs + ", " + rhs + " produced: " + DecVector3.SqrMagnitude(lhs - rhs));
//            return (decimal)DecVector3.SqrMagnitude(lhs - rhs) < 0.0 / 1.0;
        }

        public static bool operator !=(DecVector3 lhs, DecVector3 rhs) {
            return (decimal)DecVector3.SqrMagnitude(lhs - rhs) >= 0.0 / 1.0;
        }

		public static implicit operator DecVector3(Vector3 vector3) {
			return new DecVector3(vector3);
		}
        public static explicit operator Vector3(DecVector3 vector3d) {
            return new Vector3((float)vector3d.x, (float)vector3d.y, (float)vector3d.z);
        }

        public static DecVector3 Lerp(DecVector3 from, DecVector3 to, decimal t) {
            t = MathDec.Clamp01(t);
            return new DecVector3(from.x + (to.x - from.x) * t, from.y + (to.y - from.y) * t, from.z + (to.z - from.z) * t);
        }

        public static DecVector3 Slerp(DecVector3 from, DecVector3 to, decimal t) {
            Vector3 v3 = Vector3.Slerp((Vector3)from, (Vector3)to, (float)t);
            return new DecVector3(v3);
        }

        public static void OrthoNormalize(ref DecVector3 normal, ref DecVector3 tangent) {
            Vector3 v3normal = new Vector3();
            Vector3 v3tangent = new Vector3();
            v3normal = (Vector3)normal;
            v3tangent = (Vector3)tangent;
            Vector3.OrthoNormalize(ref v3normal, ref v3tangent);
            normal = new DecVector3(v3normal);
            tangent = new DecVector3(v3tangent);
        }

        public static void OrthoNormalize(ref DecVector3 normal, ref DecVector3 tangent, ref DecVector3 binormal) {
            Vector3 v3normal = new Vector3();
            Vector3 v3tangent = new Vector3();
            Vector3 v3binormal = new Vector3();
            v3normal = (Vector3)normal;
            v3tangent = (Vector3)tangent;
            v3binormal = (Vector3)binormal;
            Vector3.OrthoNormalize(ref v3normal, ref v3tangent, ref v3binormal);
            normal = new DecVector3(v3normal);
            tangent = new DecVector3(v3tangent);
            binormal = new DecVector3(v3binormal);
        }

        public static DecVector3 MoveTowards(DecVector3 current, DecVector3 target, decimal maxDistanceDelta) {
            DecVector3 vector3 = target - current;
            decimal magnitude = vector3.magnitude;
            if (magnitude <= maxDistanceDelta || magnitude == 0.0m) {
                return target;
			} else {
                return current + vector3 / magnitude * maxDistanceDelta;
			}
        }

        public static DecVector3 RotateTowards(DecVector3 current, DecVector3 target, decimal maxRadiansDelta, decimal maxMagnitudeDelta) {
            Vector3 v3 = Vector3.RotateTowards((Vector3)current, (Vector3)target, (float)maxRadiansDelta, (float)maxMagnitudeDelta);
            return new DecVector3(v3);
        }

        public static DecVector3 SmoothDamp(DecVector3 current, DecVector3 target, ref DecVector3 currentVelocity, decimal smoothTime, decimal maxSpeed) {
            decimal deltaTime = (decimal)Time.deltaTime;
            return DecVector3.SmoothDamp(current, target, ref currentVelocity, smoothTime, maxSpeed, deltaTime);
        }

        public static DecVector3 SmoothDamp(DecVector3 current, DecVector3 target, ref DecVector3 currentVelocity, decimal smoothTime) {
            decimal deltaTime = (decimal)Time.deltaTime;
            decimal maxSpeed = decimal.PositiveInfinity;
            return DecVector3.SmoothDamp(current, target, ref currentVelocity, smoothTime, maxSpeed, deltaTime);
        }

        public static DecVector3 SmoothDamp(DecVector3 current, DecVector3 target, ref DecVector3 currentVelocity, decimal smoothTime, decimal maxSpeed, decimal deltaTime) {
            smoothTime = MathDec.Max(0.0001m, smoothTime);
            decimal num1 = 2m / smoothTime;
            decimal num2 = num1 * deltaTime;
            decimal num3 = (1.0m / (1.0m + num2 + 0.479999989271164m * num2 * num2 + 0.234999999403954m * num2 * num2 * num2));
            DecVector3 vector = current - target;
            DecVector3 vector3_1 = target;
            decimal maxLength = maxSpeed * smoothTime;
            DecVector3 vector3_2 = DecVector3.ClampMagnitude(vector, maxLength);
            target = current - vector3_2;
            DecVector3 vector3_3 = (currentVelocity + num1 * vector3_2) * deltaTime;
            currentVelocity = (currentVelocity - num1 * vector3_3) * num3;
            DecVector3 vector3_4 = target + (vector3_2 + vector3_3) * num3;
            if (DecVector3.Dot(vector3_1 - current, vector3_4 - vector3_1) > 0.0) {
                vector3_4 = vector3_1;
                currentVelocity = (vector3_4 - vector3_1) / deltaTime;
            }
            return vector3_4;
        }

        public void Set(decimal new_x, decimal new_y, decimal new_z) {
            this.x = new_x;
            this.y = new_y;
            this.z = new_z;
        }

        public static DecVector3 Scale(DecVector3 a, DecVector3 b) {
            return new DecVector3(a.x * b.x, a.y * b.y, a.z * b.z);
        }

		public void Scale(Vector3 scale) {
			this.x *= scale.x;
			this.y *= scale.y;
			this.z *= scale.z;
		}
        public void Scale(DecVector3 scale) {
            this.x *= scale.x;
            this.y *= scale.y;
            this.z *= scale.z;
        }
		public void InverseScale(Vector3 scale){
			this.x /= scale.x;
			this.y /= scale.y;
			this.z /= scale.z;
		}
		public void InverseScale(DecVector3 scale){
			this.x /= scale.x;
			this.y /= scale.y;
			this.z /= scale.z;
		}

        public static DecVector3 Cross(DecVector3 lhs, DecVector3 rhs) {
            return new DecVector3(lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x);
        }

		/**
		 *	Tests if the values of the vector are approximately equal. The value of precisionBase is a multiple of the
		 *	base precision of a decimal, which is itself multiplied by the values to determine the tolerance. Any value
		 *	less than zero means "must be exactly zero".
		 */
		public static bool Approximately(DecVector3 lhs, DecVector3 rhs, decimal toleranceMultiplier){
			return MathDec.ApproximatelyRelative(lhs.x, rhs.x, toleranceMultiplier) 
					&& MathDec.ApproximatelyRelative(lhs.y, rhs.y, toleranceMultiplier)
					&& MathDec.ApproximatelyRelative(lhs.z, rhs.z, toleranceMultiplier);
		}

		public static bool ApproximatelyFixed(DecVector3 lhs, DecVector3 rhs, decimal toleranceBase){
			return MathDec.ApproximatelyFixed(lhs.x, rhs.x, toleranceBase) 
					&& MathDec.ApproximatelyFixed(lhs.y, rhs.y, toleranceBase)
					&& MathDec.ApproximatelyFixed(lhs.z, rhs.z, toleranceBase);
		}

        public override int GetHashCode() {
            return this.x.GetHashCode() ^ this.y.GetHashCode() << 2 ^ this.z.GetHashCode() >> 2;
        }

        public override bool Equals(object other) {
            if (!(other is DecVector3))
                return false;
            DecVector3 vector3d = (DecVector3)other;
            if (this.x.Equals(vector3d.x) && this.y.Equals(vector3d.y)) {
                return this.z.Equals(vector3d.z);
			} else {
                return false;
			}
        }

        public static DecVector3 Reflect(DecVector3 inDirection, DecVector3 inNormal) {
            return -2m * DecVector3.Dot(inNormal, inDirection) * inNormal + inDirection;
        }

        public static DecVector3 Normalize(DecVector3 value) {
            decimal num = DecVector3.Magnitude(value);
            if (num > 9.99999974737875E-06) {
                return value / num;
			} else {
                return DecVector3.zero;
			}
        }

        public void Normalize() {
            decimal num = DecVector3.Magnitude(this);
            if (num > 9.99999974737875E-06) {
                this = this / num;
			} else {
                this = DecVector3.zero;
			}
        }

        public override string ToString() {
			return ToString (defaultNumberFormat);
//            return "(" + this.x + ", " + this.y + ", " + this.z + ")";
        }

		public string ToString(String numberFormat) {
			return "(" + this.x.ToString(numberFormat) + ", " 
				+ this.y.ToString(numberFormat) + ", " 
				+ this.z.ToString(numberFormat) + ")";
		}

        public static decimal Dot(DecVector3 lhs, DecVector3 rhs) {
            return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
        }

		//	[Alex:] I don't know what this means. The check against an absolute value disturbs me, though.
        public static DecVector3 Project(DecVector3 vector, DecVector3 onNormal) {
            decimal num = DecVector3.Dot(onNormal, onNormal);
            if (num < 1.40129846432482E-45m) {
                return DecVector3.zero;
			} else {
                return onNormal * DecVector3.Dot(vector, onNormal) / num;
			}
        }

        public static DecVector3 Exclude(DecVector3 excludeThis, DecVector3 fromThat) {
            return fromThat - DecVector3.Project(fromThat, excludeThis);
        }

        public static decimal Angle(DecVector3 from, DecVector3 to) {
			return MathDec.Acos(MathDec.Clamp(DecVector3.Dot(from.normalized, to.normalized), -1m, 1m)) * degreesPerRadian;
        }

        public static decimal Distance(DecVector3 a, DecVector3 b) {
            DecVector3 vector3d = new DecVector3(a.x - b.x, a.y - b.y, a.z - b.z);
            return Math.Sqrt(vector3d.x * vector3d.x + vector3d.y * vector3d.y + vector3d.z * vector3d.z);

        }

        public static DecVector3 ClampMagnitude(DecVector3 vector, decimal maxLength) {
			//	Why square it? Why not use absolute value or just check for negatives?
			if (vector.sqrMagnitude > maxLength * maxLength) {
                return vector.normalized * maxLength;
			} else {
                return vector;
			}
        }

        public static decimal Magnitude(DecVector3 a) {
            return Math.Sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
        }

        public static decimal SqrMagnitude(DecVector3 a) {
            return a.x * a.x + a.y * a.y + a.z * a.z;
        }

        public static DecVector3 Min(DecVector3 lhs, DecVector3 rhs) {
            return new DecVector3(MathDec.Min(lhs.x, rhs.x), MathDec.Min(lhs.y, rhs.y), MathDec.Min(lhs.z, rhs.z));
        }

        public static DecVector3 Max(DecVector3 lhs, DecVector3 rhs) {
            return new DecVector3(MathDec.Max(lhs.x, rhs.x), MathDec.Max(lhs.y, rhs.y), MathDec.Max(lhs.z, rhs.z));
        }
		//	Non-static versions
		public DecVector3 Min(DecVector3 rhs) {
			return new DecVector3(MathDec.Min(this.x, rhs.x), MathDec.Min(this.y, rhs.y), MathDec.Min(this.z, rhs.z));
		}
		public DecVector3 Max(DecVector3 rhs) {
			return new DecVector3(MathDec.Max(this.x, rhs.x), MathDec.Max(this.y, rhs.y), MathDec.Max(this.z, rhs.z));
		}
        
    }
}

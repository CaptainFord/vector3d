// Type: UnityEngine.Vector2
// Assembly: UnityEngine, Version=0.0.0.0, Culture=neutral, PublicKeyToken=null
// Assembly location: C:\Program Files (x86)\Unity\Editor\Data\Managed\UnityEngine.dll
using System;
using System.Runtime.CompilerServices;

namespace UnityEngine {
    public struct DecVector2 {
        public const decimal kEpsilon = 1E-05m;
		public const decimal degreesPerRadian = 57.2957795130823m;
		public const string defaultNumberFormat = "G5";
        public decimal x;
        public decimal y;

        public decimal this[int index] {
            get {
                switch (index) {
                    case 0:
                        return this.x;
                    case 1:
                        return this.y;
                    default:
                        throw new IndexOutOfRangeException("Invalid DecVector2 index! (" + index + ")");
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
                    default:
						throw new IndexOutOfRangeException("Invalid DecVector2 index! (" + index + ")");
                }
            }
        }

        public DecVector2 normalized {
            get {
                DecVector2 vector2d = new DecVector2(this.x, this.y);
                vector2d.Normalize();
                return vector2d;
            }
        }

        public decimal magnitude {
            get {
                return MathDec.Sqrt(this.x * this.x + this.y * this.y);
            }
        }

        public decimal sqrMagnitude {
            get {
                return this.x * this.x + this.y * this.y;
            }
        }

        public static DecVector2 zero {
            get {
                return new DecVector2(0.0m, 0.0m);
            }
        }

        public static DecVector2 one {
            get {
                return new DecVector2(1m, 1m);
            }
        }

        public static DecVector2 up {
            get {
                return new DecVector2(0.0m, 1m);
            }
        }
		public static DecVector2 down {
			get {
				return new DecVector2(0.0m, -1m);
			}
		}

        public static DecVector2 right {
            get {
                return new DecVector2(1m, 0.0m);
            }
        }
		public static DecVector2 left {
			get {
				return new DecVector2(-1m, 0.0m);
			}
		}

        public DecVector2(decimal x, decimal y) {
            this.x = x;
            this.y = y;
        }

        public static implicit operator DecVector2(DecVector3 v) {
            return new DecVector2(v.x, v.y);
        }

        public static implicit operator DecVector3(DecVector2 v) {
            return new DecVector3(v.x, v.y, 0.0m);
        }

        public static DecVector2 operator +(DecVector2 a, DecVector2 b) {
            return new DecVector2(a.x + b.x, a.y + b.y);
        }

        public static DecVector2 operator -(DecVector2 a, DecVector2 b) {
            return new DecVector2(a.x - b.x, a.y - b.y);
        }

        public static DecVector2 operator -(DecVector2 a) {
            return new DecVector2(-a.x, -a.y);
        }

        public static DecVector2 operator *(DecVector2 a, decimal d) {
            return new DecVector2(a.x * d, a.y * d);
        }

        public static DecVector2 operator *(float d, DecVector2 a) {
            return new DecVector2(a.x * d, a.y * d);
        }

        public static DecVector2 operator /(DecVector2 a, decimal d) {
            return new DecVector2(a.x / d, a.y / d);
        }

        public static bool operator ==(DecVector2 lhs, DecVector2 rhs) {
            return DecVector2.SqrMagnitude(lhs - rhs) < 0.0m / 1.0m;
        }

        public static bool operator !=(DecVector2 lhs, DecVector2 rhs) {
            return (decimal)DecVector2.SqrMagnitude(lhs - rhs) >= 0.0m / 1.0m;
        }

        public void Set(decimal new_x, decimal new_y) {
            this.x = new_x;
            this.y = new_y;
        }

        public static DecVector2 Lerp(DecVector2 from, DecVector2 to, decimal t) {
            t = MathDec.Clamp01(t);
            return new DecVector2(from.x + (to.x - from.x) * t, from.y + (to.y - from.y) * t);
        }

        public static DecVector2 MoveTowards(DecVector2 current, DecVector2 target, decimal maxDistanceDelta) {
            DecVector2 vector2 = target - current;
            decimal magnitude = vector2.magnitude;
            if (magnitude <= maxDistanceDelta || magnitude == 0.0m)
                return target;
            else
                return current + vector2 / magnitude * maxDistanceDelta;
        }

        public static DecVector2 Scale(DecVector2 a, DecVector2 b) {
            return new DecVector2(a.x * b.x, a.y * b.y);
        }

        public void Scale(DecVector2 scale) {
            this.x *= scale.x;
            this.y *= scale.y;
        }

        public void Normalize() {
            decimal magnitude = this.magnitude;
            if (magnitude > 9.99999974737875E-06) {
                this = this / magnitude;
			} else {
                this = DecVector2.zero;
			}
        }

        public override string ToString() {
			return ToString (defaultNumberFormat);
		}

        public string ToString(string format) {
			return "(" + x.ToString(format) + "," + y.ToString(format) + ")";
        }

        public override int GetHashCode() {
            return this.x.GetHashCode() ^ this.y.GetHashCode() << 2;
        }

        public override bool Equals(object other) {
            if (!(other is DecVector2))
                return false;
            DecVector2 vector2d = (DecVector2)other;
            if (this.x.Equals(vector2d.x))
                return this.y.Equals(vector2d.y);
            else
                return false;
        }

        public static decimal Dot(DecVector2 lhs, DecVector2 rhs) {
            return lhs.x * rhs.x + lhs.y * rhs.y;
        }

        public static decimal Angle(DecVector2 from, DecVector2 to) {
			return MathDec.Acos(MathDec.Clamp(DecVector2.Dot(from.normalized, to.normalized), -1m, 1m)) * degreesPerRadian;
        }

        public static decimal Distance(DecVector2 a, DecVector2 b) {
            return (a - b).magnitude;
        }

        public static DecVector2 ClampMagnitude(DecVector2 vector, decimal maxLength) {
            if (vector.sqrMagnitude > maxLength * maxLength) {
                return vector.normalized * maxLength;
			} else {
                return vector;
			}
        }

        public static decimal SqrMagnitude(DecVector2 a) {
            return (a.x * a.x + a.y * a.y);
        }

        public decimal SqrMagnitude() {
            return (this.x * this.x + this.y * this.y);
        }


        public static DecVector2 Min(DecVector2 lhs, DecVector2 rhs) {
            return new DecVector2(MathDec.Min(lhs.x, rhs.x), MathDec.Min(lhs.y, rhs.y));
        }

        public static DecVector2 Max(DecVector2 lhs, DecVector2 rhs) {
            return new DecVector2(MathDec.Max(lhs.x, rhs.x), MathDec.Max(lhs.y, rhs.y));
        }
		//	Why do it as static? I guess that kind of makes sense. But it also 
		//	makes for longer code.
		public DecVector2 Min(DecVector2 rhs) {
			return new DecVector2(MathDec.Min(x, rhs.x), MathDec.Min(y, rhs.y));
		}
		public DecVector2 Max(DecVector2 rhs) {
			return new DecVector2(MathDec.Max(x, rhs.x), MathDec.Max(y, rhs.y));
		}
    }
}

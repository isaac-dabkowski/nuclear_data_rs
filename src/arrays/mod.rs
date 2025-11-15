mod izaw;
mod jxs;
mod nxs;
mod xxs;

pub use izaw::IzawArray;
pub use jxs::JxsArray;
pub use nxs::NxsArray;
pub use xxs::XxsArray;

pub struct Arrays<'a> {
    pub nxs: &'a NxsArray,
    pub jxs: &'a JxsArray,
    pub xxs: &'a XxsArray,
}

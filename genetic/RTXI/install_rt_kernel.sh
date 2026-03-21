#! /bin/bash
set -eu

if ! id | grep -q root; then
	echo "Must run script as root; try again with sudo ./install_rt_kernel.sh."
	exit
fi

# Export environment variables
echo  "-----> Setting up variables."
export linux_version=4.19.177
export xenomai_version=3.1.3
export xenomai_root=/opt/xenomai-$xenomai_version
export build_root=/opt/build
export opt=/opt
export ipipe_patch_digit=17
export ipipe_cip_str=-cip44
export linux_tree=/opt/linux-cip-${linux_version}${ipipe_cip_str}

echo  "-----> Preparando entorno."
apt-get update && apt-get install -y autoconf automake libtool libelf-dev libssl-dev flex bison

# Descargas
cd $opt
if [ ! -f "linux-cip-${linux_version}${ipipe_cip_str}.tar.gz" ]; then
    wget --no-check-certificate https://git.kernel.org/pub/scm/linux/kernel/git/cip/linux-cip.git/snapshot/linux-cip-${linux_version}${ipipe_cip_str}.tar.gz
fi

if [ ! -f "xenomai-$xenomai_version.tar.bz2" ]; then
    wget --no-check-certificate https://source.denx.de/Xenomai/xenomai/-/archive/v$xenomai_version/xenomai-v$xenomai_version.tar.bz2 -O xenomai-$xenomai_version.tar.bz2
fi

if [ ! -f "ipipe-core-${linux_version}${ipipe_cip_str}-x86-${ipipe_patch_digit}.patch" ]; then
    wget --no-check-certificate -O ipipe-core-${linux_version}${ipipe_cip_str}-x86-${ipipe_patch_digit}.patch "https://ftp.denx.de/pub/xenomai/ipipe/v${linux_version:0:1}.x/x86/ipipe-core-${linux_version}${ipipe_cip_str}-x86-${ipipe_patch_digit}.patch"
fi

if [ ! -d "$linux_tree" ]; then
    echo "-----> Extrayendo kernel..."
    tar xf linux-cip-${linux_version}${ipipe_cip_str}.tar.gz
fi

if [ ! -d "xenomai-$xenomai_version" ]; then
    tar xf xenomai-$xenomai_version.tar.bz2
    mv xenomai-v$xenomai_version xenomai-$xenomai_version
fi

# Parchear kernel (con control para no repetir)
echo  "-----> Patching kernel."
cd $linux_tree
if [ ! -f .xenomai_prepared ]; then
    $xenomai_root/scripts/prepare-kernel.sh \
        --arch=x86 \
        --ipipe=$opt/ipipe-core-${linux_version}${ipipe_cip_str}-x86-${ipipe_patch_digit}.patch \
        --linux=$linux_tree \
        --verbose || true
    touch .xenomai_prepared
else
    echo "Kernel ya parcheado. Omitiendo."
fi

# Configuración
echo "-----> Configurando kernel (Xenomai Cobalt + Drivers)..."
if [ ! -f .config ]; then
    cp /boot/config-$(uname -r) .config
fi

# Aplicamos los flags de Xenomai
./scripts/config --enable CONFIG_IPIPE
./scripts/config --enable CONFIG_XENOMAI
./scripts/config --enable CONFIG_XENO_OPT_COBALT

# Soporte de entrada (Mouse/Trackpad)
./scripts/config --enable CONFIG_INPUT_MOUSEDEV
./scripts/config --enable CONFIG_INPUT_EVDEV
./scripts/config --enable CONFIG_INPUT_MOUSE
./scripts/config --enable CONFIG_MOUSE_PS2
./scripts/config --enable CONFIG_MOUSE_PS2_ELANTECH
./scripts/config --enable CONFIG_MOUSE_PS2_SYNAPTICS
./scripts/config --enable CONFIG_I2C_HID
./scripts/config --enable CONFIG_I2C_DESIGNWARE_PLATFORM
./scripts/config --enable CONFIG_I2C_DESIGNWARE_PCI
./scripts/config --enable CONFIG_HID_MULTITOUCH
./scripts/config --enable CONFIG_HID_GENERIC
./scripts/config --enable CONFIG_USB_HID

# Drivers de disco y sistema (Built-in)
./scripts/config --enable CONFIG_SATA_AHCI
./scripts/config --enable CONFIG_BLK_DEV_SD
./scripts/config --enable CONFIG_BLK_DEV_NVME
./scripts/config --enable CONFIG_EXT4_FS

# --- DRIVERS ESPECIFICOS (WIFI, LAN, GPU) COMPILADOS ESTÁTICAMENTE ---
./scripts/config --enable CONFIG_FW_LOADER

# Subredes y WiFi base
./scripts/config --enable CONFIG_NETDEVICES
./scripts/config --enable CONFIG_ETHERNET
./scripts/config --enable CONFIG_WLAN
./scripts/config --enable CONFIG_CFG80211
./scripts/config --enable CONFIG_MAC80211

# LAN: Intel e1000e
./scripts/config --enable CONFIG_NET_VENDOR_INTEL
./scripts/config --enable CONFIG_E1000E

# Wi-Fi: Intel iwlwifi 8260
./scripts/config --enable CONFIG_WLAN_VENDOR_INTEL
./scripts/config --enable CONFIG_IWLWIFI
./scripts/config --enable CONFIG_IWLMVM

# GPU: Intel i915 (Displays extra VGA/HDMI) y dependencias base
./scripts/config --enable CONFIG_DRM
./scripts/config --enable CONFIG_I2C_ALGOBIT
./scripts/config --enable CONFIG_DRM_I915
./scripts/config --enable CONFIG_DRM_I915_KMS
# ---------------------------------------------------------------------

# Desactivar Virtualización/Paravirtualización (Conflictos I-pipe)
./scripts/config --disable CONFIG_PARAVIRT
./scripts/config --disable CONFIG_PARAVIRT_GUEST
./scripts/config --disable CONFIG_PARAVIRT_SPINLOCKS
./scripts/config --disable CONFIG_XEN
./scripts/config --disable CONFIG_XEN_PVH
./scripts/config --disable CONFIG_KVM
./scripts/config --disable CONFIG_HYPERVISOR_GUEST
./scripts/config --disable CONFIG_HYPERV

# Neutralizar ACPI Processor, CPU Freq y CPU Idle
./scripts/config --disable CONFIG_CPU_FREQ
./scripts/config --disable CONFIG_CPU_IDLE
./scripts/config --disable CONFIG_ACPI_PROCESSOR
./scripts/config --disable CONFIG_ACPI_PROCESSOR_IDLE
./scripts/config --disable CONFIG_ACPI_CPU_FREQ
./scripts/config --disable CONFIG_INTEL_IDLE

# Desactivar Lockdep y Tracing (Causa de arch_mangle_irq_bits)
./scripts/config --disable CONFIG_PROVE_LOCKING
./scripts/config --disable CONFIG_LOCKDEP
./scripts/config --disable CONFIG_DEBUG_SPINLOCK
./scripts/config --disable CONFIG_DEBUG_MUTEXES
./scripts/config --disable CONFIG_DEBUG_LOCK_ALLOC
./scripts/config --disable CONFIG_IPIPE_TRACE
./scripts/config --disable CONFIG_TRACING
./scripts/config --disable CONFIG_FTRACE

# Desactivar firmas y debug genérico
./scripts/config --disable MODULE_SIG
./scripts/config --disable SYSTEM_TRUSTED_KEYRING
./scripts/config --disable DEBUG_INFO
./scripts/config --set-str SYSTEM_TRUSTED_KEYS ""

yes "" | make olddefconfig

# Compilación (Incremental garantizada por make al no hacer clean)
echo  "-----> Compiling kernel (Incremental)."
find tools/ -type f -name "Makefile*" -exec sed -i 's/-Werror//g' {} +
make -j$(nproc) bindeb-pkg LOCALVERSION=-xenomai-$xenomai_version KCFLAGS="-Wno-error=implicit-function-declaration -Wno-error" HOSTCFLAGS="-Wno-error"

# Instalación
echo  "-----> Installing compiled kernel"
cd $opt
dpkg -i linux-image-*xenomai-$xenomai_version*.deb
dpkg -i linux-headers-*xenomai-$xenomai_version*.deb

# Grub
echo  "-----> Configurando GRUB..."
sed -i 's/GRUB_CMDLINE_LINUX_DEFAULT=.*/GRUB_CMDLINE_LINUX_DEFAULT="quiet splash intel_idle.max_cstate=0 processor.max_cstate=1"/g' /etc/default/grub || true
update-grub

# Librerías de usuario
echo  "-----> Installing user libraries."
cd $xenomai_root
if [ ! -f configure ]; then
    ./scripts/bootstrap
fi

if [ ! -d $build_root ]; then
    mkdir -p $build_root
    cd $build_root
    $xenomai_root/configure --with-core=cobalt --enable-pshared --enable-smp --enable-dlopen-libs
else
    cd $build_root
fi

make -sj$(nproc)
make install

# Analogy y Permisos
cp -f /usr/xenomai/sbin/analogy_config /usr/sbin/
grep -q xenomai /etc/group || groupadd xenomai
usermod -a -G xenomai "$SUDO_USER"

echo  "-----> Kernel patch complete. REINICIA."
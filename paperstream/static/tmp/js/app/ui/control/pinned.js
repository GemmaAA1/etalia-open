define(['jquery', 'app/util/utils'], function ($, Utils) {

    var Pinned = function (options) {
        this.config = $.extend({
            debug: false,
            element: '#toggle-pinned'
        }, options);

        this.$element = $(this.config.element);
    };

    Pinned.prototype.init = function() {
        var that = this,
            $body = $('body');

        this.$element.on('click', function(e) {
            Utils.toggleClass(that.$element, 'active');

            $body.trigger('etalia.control.pinned.change');

            e.preventDefault();
            return false;
        });

        return this;
    };

    Pinned.prototype.setValue = function(value) {
        Utils.toggleClass(this.$element, 'active', value);

        return this;
    };

    Pinned.prototype.getValue = function() {
        return this.$element.hasClass('active');
    };

    return Pinned;
});

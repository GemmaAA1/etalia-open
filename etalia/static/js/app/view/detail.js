define([
    'app',
    'text!app/templates/detail.hbs',
    'app/model/detail'
], function (App, template) {

    var defaults = {};

    return App.View.Detail = App.Backbone.View.extend({
        tagName: 'div',
        id: 'detail',

        template: App.Handlebars.compile(template),

        events: {
            "click #detail-close": "onCloseClick",
            "click #detail-prev": "onPrevClick",
            "click #detail-next": "onNextClick",
            "click": 'onClick'
        },

        initialize: function (options) {
            App.defaults(options, defaults);

            if (!this.model instanceof App.Model.Detail) {
                throw 'Unexpected model type';
            }
        },

        render: function () {
            var prev = this.model.get('prev'),
                next = this.model.get('next');

            App.$('div[data-detail-placeholder]').replaceWith(this.$el.html(this.template({
                prev: prev ? prev.attributes : null,
                next: next ? next.attributes : null
            })));

            var that = this,
                view = this.model.get('view').render();
            view.$el.appendTo(this.$('.document .wrapper'));

            this.listenToOnce(view, 'close', function() {
                that.close();
            });

            return this.open();
        },

        onPrevClick: function() {
            var prev = this.model.get('prev');
            if (prev) {
                this.trigger('detail:prev', prev, this);
            }
        },

        onNextClick: function() {
            var next = this.model.get('next');
            if (next) {
                this.trigger('detail:next', next, this);
            }
        },

        onClick: function(e) {
            if (App.$(e.target).attr('id') === 'detail-document') {
                this.close();
            }
        },

        onCloseClick: function() {
            this.close();
        },

        open: function() {
            App.$('body').addClass('detail-opened');

            return this;
        },

        close: function() {
            App.$('body').removeClass('detail-opened');

            this.model.get('view').remove();

            return this;
        }
    });
});
